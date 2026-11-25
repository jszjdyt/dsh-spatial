#!/usr/bin/env Rscript
# dsh-spatial — RCTD（spacexr）官方实现脚本
#
# 调用：Rscript --vanilla lib/r_rctd.R --in <bundle> --out <out> --params '<json>'
#
# bundle（由 lib/sp_rbridge.py 的 build_rctd_bundle 写出）：
#   ref_counts.tsv.gz | ref_counts.mtx          参考 counts（gene × cell，原始整数 counts）
#   ref_types.tsv                               参考细胞类型（cell, cell_type 两列）
#   spot_counts.tsv.gz | spot_counts.mtx        空间 counts（gene × spot）
#   coords.tsv                                  坐标（spot, x, y）
#   genes.txt                                   mtx 时的行名（gene，顺序与 mtx 行一致）
#   params.json                                 mode / max_cores / min_genes / cell_type_col 等
#
# 输出：
#   <out>/result.json                           method/mode/cell_types/n_spots/n_cells/warnings 等
#   <out>/tables/weights.csv                    spot × cell_type 比例（行和=1，spacexr::normalize_weights 归一化；
#                                               full 模式的原始权重只保证 ≤1，见 spacexr 官方教程）
#   <out>/tables/weights_raw.csv                未归一化的原始权重（可审计）
#   <out>/tables/results_df.csv                 doublet/multi 的 spot 分类；full 模式该表为空（spacexr 全混合模式
#                                               不产生 spot 分类，属正常行为，不伪造分类）
#   <out>/error.json                            出错时的 {code, message}
#
# 退出码：0 成功；3 缺前置（数据不满足，转 prerequisite）；4 运行期程序错误（上层不得据此改换方法）；
#         5 缺依赖包（转 unsupported）；其它 失败。
#
# 真实实现声明：本脚本的全部估计来自 spacexr::Reference / spacexr::SpatialRNA /
# spacexr::create.RCTD / spacexr::run.RCTD（RCTD 论文 Cable et al. 2022 Nat Biotechnol），
# 只负责拆包/打包与导出；不做近似重算、不冒充其它方法、不静默改变 doublet_mode。
# 权重归一化用 spacexr 自带的 normalize_weights()（spacexr 教程里 full 模式的标准后处理），
# 原始权重另存 weights_raw.csv 以便复核。

suppressPackageStartupMessages({
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    cat("MISSING_PACKAGE:spacexr\n", file = stderr())
    cat("R 包 spacexr（RCTD）未安装：无法调用 spacexr::create.RCTD / run.RCTD 官方实现\n", file = stderr())
    quit(save = "no", status = 5)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    cat("MISSING_PACKAGE:Matrix\n", file = stderr())
    cat("R 包 Matrix 未安装：无法读取 mtx 计数矩阵\n", file = stderr())
    quit(save = "no", status = 5)
  }
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(key) {
  idx <- match(key, args)
  if (is.na(idx) || idx >= length(args)) return(NULL)
  args[idx + 1]
}

in_dir <- getarg("--in")
out_dir <- getarg("--out")
params_raw <- getarg("--params")

if (is.null(in_dir) || is.null(out_dir)) {
  cat("R 脚本缺少 --in / --out 参数\n", file = stderr())
  quit(save = "no", status = 3)
}
if (!dir.exists(in_dir)) {
  cat(paste0("bundle 目录不存在：", in_dir, "\n"), file = stderr())
  quit(save = "no", status = 3)
}
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 极简 JSON（读写都不依赖 jsonlite：spacexr 的依赖里没有它）
# ---------------------------------------------------------------------------

json_escape <- function(s) {
  s <- enc2utf8(as.character(s))
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s
}

json_emit <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    nm <- names(x)
    if (!is.null(nm) && length(nm) == length(x) && all(nzchar(nm))) {
      parts <- vapply(seq_along(x), function(i) {
        paste0("\"", json_escape(nm[i]), "\":", json_emit(x[[i]]))
      }, character(1))
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    parts <- vapply(seq_along(x), function(i) json_emit(x[[i]]), character(1))
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (length(x) == 1L) {
    if (is.logical(x) && is.na(x)) return("null")
    if (is.logical(x)) return(if (isTRUE(x)) "true" else "false")
    if (is.numeric(x)) {
      if (!is.finite(x)) return("null")
      return(format(x, scientific = FALSE, trim = TRUE, digits = 12))
    }
    if (is.na(x)) return("null")
    return(paste0("\"", json_escape(x), "\""))
  }
  parts <- vapply(seq_along(x), function(i) json_emit(x[[i]]), character(1))
  paste0("[", paste(parts, collapse = ","), "]")
}

write_json <- function(obj, path) {
  writeLines(json_emit(obj), path, useBytes = TRUE)
}

json_unescape <- function(s) {
  s <- gsub("\\\"", "\"", s, fixed = TRUE)
  s <- gsub("\\\\", "\\", s, fixed = TRUE)
  s <- gsub("\\n", "\n", s, fixed = TRUE)
  s <- gsub("\\t", "\t", s, fixed = TRUE)
  s
}

#: 解析扁平 JSON 对象（值只允许字符串/数字/true/false/null）；解析不了的键直接忽略。
parse_params <- function(raw) {
  out <- list()
  if (is.null(raw) || !nzchar(trimws(raw))) return(out)
  txt <- trimws(raw)
  if (substr(txt, 1L, 1L) != "{" || substr(txt, nchar(txt), nchar(txt)) != "}") return(out)
  body <- substr(txt, 2L, nchar(txt) - 1L)
  chars <- strsplit(body, "", fixed = TRUE)[[1]]
  parts <- character(0)
  buf <- ""
  in_str <- FALSE
  escaped <- FALSE
  depth <- 0L
  for (ch in chars) {
    if (in_str) {
      buf <- paste0(buf, ch)
      if (escaped) {
        escaped <- FALSE
      } else if (ch == "\\") {
        escaped <- TRUE
      } else if (ch == "\"") {
        in_str <- FALSE
      }
      next
    }
    if (ch == "\"") {
      in_str <- TRUE
      buf <- paste0(buf, ch)
      next
    }
    if (ch == "{" || ch == "[") depth <- depth + 1L
    if (ch == "}" || ch == "]") depth <- depth - 1L
    if (ch == "," && depth == 0L) {
      parts <- c(parts, buf)
      buf <- ""
      next
    }
    buf <- paste0(buf, ch)
  }
  if (nzchar(trimws(buf))) parts <- c(parts, buf)
  for (part in parts) {
    m <- regmatches(part, regexec("^\\s*\"((?:[^\"\\\\]|\\\\.)*)\"\\s*:\\s*(.*)$", part, perl = TRUE))[[1]]
    if (length(m) != 3L) next
    key <- json_unescape(m[2])
    value <- trimws(m[3])
    if (substr(value, 1L, 1L) == "\"" && nchar(value) >= 2L) {
      out[[key]] <- json_unescape(substr(value, 2L, nchar(value) - 1L))
    } else if (value == "true") {
      out[[key]] <- TRUE
    } else if (value == "false") {
      out[[key]] <- FALSE
    } else if (value == "null") {
      out[[key]] <- NULL
    } else {
      num <- suppressWarnings(as.numeric(value))
      if (!is.na(num)) out[[key]] <- num
    }
  }
  out
}

# ---------------------------------------------------------------------------
# 错误写回 + 退出码（与 sp_rbridge.py 的映射一致）
# ---------------------------------------------------------------------------

write_error <- function(code, message) {
  try(
    write_json(list(code = code, message = as.character(message)),
               file.path(out_dir, "error.json")),
    silent = TRUE
  )
  # 统一写 UTF-8 字节：上层 Python 用 utf-8 读 error.json / stderr 才不会乱码
  writeLines(enc2utf8(as.character(message)), con = stderr(), useBytes = TRUE)
}

#: 缺前置（数据不满足）：上层转 prerequisite
fail_prereq <- function(message) {
  write_error("prerequisite", message)
  quit(save = "no", status = 3)
}

#: 运行期程序错误：上层**不得**据此改换反卷积方法
fail_runtime <- function(message) {
  write_error("runtime", message)
  quit(save = "no", status = 4)
}

# ---------------------------------------------------------------------------
# 参数
# ---------------------------------------------------------------------------

params <- parse_params(params_raw)
pick <- function(key, default) {
  value <- params[[key]]
  if (is.null(value) || (length(value) == 1L && is.na(value))) return(default)
  value
}
pick_num <- function(key, default) {
  value <- suppressWarnings(as.numeric(pick(key, default)))
  if (length(value) != 1L || is.na(value)) default else value
}

mode <- as.character(pick("mode", "full"))
if (!(mode %in% c("full", "doublet", "multi"))) {
  fail_prereq(sprintf("params.mode=%s 不合法：只能是 full / doublet / multi", mode))
}
min_genes <- as.integer(pick_num("min_genes", 50))
if (is.na(min_genes) || min_genes < 1L) min_genes <- 50L
cell_type_col <- as.character(pick("cell_type_col", ""))
max_cores <- as.integer(pick_num("max_cores", 1))
if (is.na(max_cores)) max_cores <- 1L
max_cores <- max(1L, min(max_cores, 4L))  # 契约默认 1；最多 4，避免把机器打满
cell_min <- as.integer(pick_num("cell_min", 25))
if (is.na(cell_min) || cell_min < 3L) cell_min <- 25L

WARNINGS <- character(0)
add_warn <- function(msg) {
  msg <- as.character(msg)
  if (!length(msg)) return(invisible(NULL))
  msg <- msg[1]
  if (length(WARNINGS) < 50L && !(msg %in% WARNINGS)) WARNINGS <<- c(WARNINGS, msg)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# 读 bundle
# ---------------------------------------------------------------------------

#: read.delim 会按魔数自动解压 .gz/.bz2/.xz（见 ?read.table 的 file 参数），
#: 不要再套 gzfile()：未显式打开的连接在这里会报 "invalid connection"。
read_delim_file <- function(path) {
  read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("NA", "NaN", ""))
}

#: 读 gene × sample 计数矩阵：优先 mtx（大矩阵省内存/时间），其次 tsv.gz/tsv。
#: mtx 不带行/列名：行名取自 genes.txt，列名由调用方按 ref_types/coords 的顺序传入。
read_counts_matrix <- function(dir, base, expected_cols = NULL) {
  mtx <- file.path(dir, paste0(base, ".mtx"))
  genes_path <- file.path(dir, "genes.txt")
  if (file.exists(mtx)) {
    if (!file.exists(genes_path)) {
      fail_prereq(sprintf("bundle 缺 genes.txt：%s.mtx 需要行名文件", base))
    }
    mat <- tryCatch(as.matrix(Matrix::readMM(mtx)), error = function(e) NULL)
    if (is.null(mat)) fail_prereq(sprintf("%s.mtx 读取失败（Matrix::readMM）", base))
    genes <- readLines(genes_path, warn = FALSE)
    if (length(genes) != nrow(mat)) {
      fail_prereq(sprintf("genes.txt 行数（%d）与 %s.mtx 行数（%d）不一致", length(genes), base, nrow(mat)))
    }
    rownames(mat) <- genes
    if (!is.null(expected_cols) && length(expected_cols) == ncol(mat)) {
      colnames(mat) <- expected_cols
    }
    return(mat)
  }
  for (suffix in c(".tsv.gz", ".tsv", ".txt.gz", ".txt")) {
    path <- file.path(dir, paste0(base, suffix))
    if (file.exists(path)) {
      mat <- tryCatch(as.matrix(read_delim_file(path)), error = function(e) e)
      if (inherits(mat, "error")) {
        fail_prereq(sprintf("%s%s 读取失败：%s", base, suffix, conditionMessage(mat)))
      }
      return(mat)
    }
  }
  fail_prereq(sprintf("bundle 里找不到 %s.mtx / %s.tsv.gz / %s.tsv", base, base, base))
}

need_file <- function(name, hint) {
  path <- file.path(in_dir, name)
  if (!file.exists(path)) fail_prereq(sprintf("bundle 缺 %s（%s）", name, hint))
  path
}

ref_types <- tryCatch(read.delim(need_file("ref_types.tsv", "cell, cell_type 两列"),
                                 stringsAsFactors = FALSE, check.names = FALSE),
                      error = function(e) NULL)
coords_df <- tryCatch(read.delim(need_file("coords.tsv", "spot, x, y 三列"),
                                 stringsAsFactors = FALSE, check.names = FALSE),
                      error = function(e) NULL)
if (is.null(ref_types) || !all(c("cell", "cell_type") %in% colnames(ref_types))) {
  fail_prereq("ref_types.tsv 读取失败或缺少 cell / cell_type 列")
}

ref_mat <- read_counts_matrix(in_dir, "ref_counts", expected_cols = ref_types$cell)
spot_cols <- if (!is.null(coords_df) && nrow(coords_df) > 0L && "spot" %in% colnames(coords_df)) {
  coords_df$spot
} else {
  NULL
}
spot_mat <- read_counts_matrix(in_dir, "spot_counts", expected_cols = spot_cols)
if (is.null(colnames(spot_mat))) {
  # mtx + 没有 coords.tsv 行：spot 名字无处可取，用占位名（下游只用名字做对齐）
  colnames(spot_mat) <- paste0("spot", seq_len(ncol(spot_mat)))
  add_warn("coords.tsv 里没有 spot 行名：spot 用占位名 spot1..spotN（只影响标签，不影响权重）")
}

# ---------------------------------------------------------------------------
# 数据前置检查（退出码 3）
# ---------------------------------------------------------------------------

check_counts <- function(mat, label) {
  vals <- as.numeric(mat)
  if (!length(vals)) fail_prereq(sprintf("%s 是空矩阵", label))
  if (any(!is.finite(vals))) fail_prereq(sprintf("%s 含 NA/Inf：RCTD 需要原始 counts", label))
  if (any(vals < 0)) fail_prereq(sprintf("%s 含负值：RCTD 需要原始非负 counts", label))
  bad <- abs(vals - round(vals)) > 1e-6
  if (any(bad)) {
    fail_prereq(sprintf(
      "%s 有 %.1f%% 的值不是整数（max=%.4g）：RCTD 需要原始 counts，不接受 log/标准化后的矩阵",
      label, 100 * mean(bad), max(vals)))
  }
}

check_counts(ref_mat, "参考矩阵 ref_counts")
check_counts(spot_mat, "空间矩阵 spot_counts")

# 细胞名对齐（ref_types.tsv → ref_counts 的列）
if (!setequal(colnames(ref_mat), ref_types$cell)) {
  matched <- length(intersect(colnames(ref_mat), ref_types$cell))
  fail_prereq(sprintf("ref_counts 的列名与 ref_types.tsv 的 cell 名不一致（只匹配 %d/%d）",
                      matched, ncol(ref_mat)))
}
ref_types <- ref_types[match(colnames(ref_mat), ref_types$cell), , drop = FALSE]

# spot 名对齐（coords.tsv → spot_counts 的列）
if (is.null(coords_df) || nrow(coords_df) == 0L ||
    !all(c("spot", "x", "y") %in% colnames(coords_df))) {
  add_warn("coords.tsv 为空：用 (0,0) 占位坐标（RCTD full 模式的权重不依赖坐标）")
  coords_used <- data.frame(spot = colnames(spot_mat), x = 0, y = 0,
                            stringsAsFactors = FALSE)
  coords_kind <- "fake"
} else {
  if (!setequal(colnames(spot_mat), coords_df$spot)) {
    matched <- length(intersect(colnames(spot_mat), coords_df$spot))
    fail_prereq(sprintf("spot_counts 的列名与 coords.tsv 的 spot 名不一致（只匹配 %d/%d）",
                        matched, ncol(spot_mat)))
  }
  coords_df <- coords_df[match(colnames(spot_mat), coords_df$spot), , drop = FALSE]
  coords_used <- data.frame(spot = coords_df$spot,
                            x = suppressWarnings(as.numeric(coords_df$x)),
                            y = suppressWarnings(as.numeric(coords_df$y)),
                            stringsAsFactors = FALSE)
  if (any(!is.finite(coords_used$x)) || any(!is.finite(coords_used$y))) {
    add_warn("coords.tsv 里有非数值坐标：用 (0,0) 占位")
    coords_used$x[!is.finite(coords_used$x)] <- 0
    coords_used$y[!is.finite(coords_used$y)] <- 0
    coords_kind <- "fake"
  } else {
    coords_kind <- "real"
  }
}

# 共同基因（RCTD 要求参考与空间矩阵同一套基因；Python 侧已取交集，这里再兜底校验）
common <- intersect(rownames(ref_mat), rownames(spot_mat))
if (length(common) < min_genes) {
  fail_prereq(sprintf(
    "参考与空间数据的共同基因只有 %d 个 < min_genes=%d：请确认两者基因命名一致（symbol/Ensembl），或降低 min_genes",
    length(common), min_genes))
}
ref_mat <- ref_mat[common, , drop = FALSE]
spot_mat <- spot_mat[common, , drop = FALSE]

# 细胞类型
types <- ref_types$cell_type
cell_types <- sort(unique(as.character(types)))
if (length(cell_types) < 2L) {
  fail_prereq(sprintf("参考里只有 %d 种细胞类型：RCTD 至少需要 2 种（列 %s）",
                      length(cell_types), ifelse(nzchar(cell_type_col), cell_type_col, "cell_type")))
}
type_counts <- table(factor(types, levels = cell_types))
if (min(type_counts) < cell_min) {
  few <- names(type_counts)[which.min(type_counts)]
  fail_prereq(sprintf(
    "参考里细胞类型「%s」只有 %d 个细胞 < %d（RCTD 的 CELL_MIN_INSTANCE）：请补足该类型细胞，或下调 cell_min 参数",
    few, as.integer(min(type_counts)), cell_min))
}
ct <- factor(as.character(types), levels = cell_types)
names(ct) <- colnames(ref_mat)

# 空间侧 UMI 太少时 create.RCTD 会丢掉 spot：提前提示（不阻断）
n_umi <- colSums(spot_mat)
if (any(n_umi <= 0)) {
  add_warn(sprintf("有 %d 个 spot 的 counts 全为 0：RCTD 会丢弃这些 spot", sum(n_umi <= 0)))
}

# ---------------------------------------------------------------------------
# RCTD（真实调用 spacexr）
# ---------------------------------------------------------------------------

t_start <- Sys.time()
run_one <- function(expr) {
  withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      add_warn(conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

#: spacexr 里属于"数据不够/不匹配"的错误走 prerequisite（退出码 3），
#: 真正实现不了的错误才走 runtime（退出码 4，上层不得据此改换方法）。
fail_rctd <- function(prefix, msg) {
  msg <- as.character(msg)
  if (grepl("fewer than 10 .*differentially expressed genes", msg)) {
    fail_prereq(paste0(
      prefix,
      "：参考里细胞类型之间的差异基因不足 10 个（spacexr 的 get_de_genes 要求 ≥10 个用于回归/平台效应拟合）：",
      "请换细胞类型区分度更高的参考（例如更细的 cell_type），或检查 cell_type 列是否把不同类型混在一起；",
      "spacexr 原始信息：", msg))
  }
  if (grepl("At least 90% of genes do not match", msg)) {
    fail_prereq(paste0(
      prefix,
      "：参考与空间数据的基因名匹配率不足 10%（命名体系不一致，如 symbol vs Ensembl）：请统一注释后重试；",
      "spacexr 原始信息：", msg))
  }
  if (grepl("need a minimum of .*cells for each cell type", msg)) {
    fail_prereq(paste0(
      prefix, "：参考里每个细胞类型的细胞数不足（CELL_MIN_INSTANCE）：请补足细胞或下调 cell_min；",
      "spacexr 原始信息：", msg))
  }
  if (grepl("no cell types passed weight threshold", msg)) {
    fail_prereq(paste0(
      prefix, "：没有细胞类型的权重超过初筛阈值（通常是 spot counts 太少或参考谱完全不匹配）：",
      "请检查空间矩阵是否原始 counts、参考与数据是否同一组织；spacexr 原始信息：", msg))
  }
  fail_runtime(paste0(prefix, "：", msg))
}

ref_obj <- run_one(spacexr::Reference(counts = ref_mat, cell_types = ct,
                                      nUMI = colSums(ref_mat)))
if (inherits(ref_obj, "error")) {
  fail_rctd("spacexr::Reference 构造失败", conditionMessage(ref_obj))
}

sp_obj <- run_one(spacexr::SpatialRNA(
  coords = data.frame(x = coords_used$x, y = coords_used$y, row.names = coords_used$spot),
  counts = spot_mat,
  nUMI = n_umi))
if (inherits(sp_obj, "error")) {
  fail_rctd("spacexr::SpatialRNA 构造失败", conditionMessage(sp_obj))
}

# create.RCTD 的可选参数：只在调用方显式给了才覆盖 spacexr 默认值，避免"偷偷改方法"
create_args <- list(
  spatialRNA = sp_obj,
  reference = ref_obj,
  max_cores = max_cores,
  test_mode = FALSE,
  CELL_MIN_INSTANCE = cell_min
)
optional <- list(
  UMI_min = params$umi_min,
  UMI_max = params$umi_max,
  counts_MIN = params$counts_min,
  gene_cutoff = params$gene_cutoff,
  fc_cutoff = params$fc_cutoff,
  gene_cutoff_reg = params$gene_cutoff_reg,
  fc_cutoff_reg = params$fc_cutoff_reg,
  MAX_MULTI_TYPES = params$max_multi_types,
  CONFIDENCE_THRESHOLD = params$confidence_threshold,
  DOUBLET_THRESHOLD = params$doublet_threshold
)
for (name in names(optional)) {
  if (!is.null(optional[[name]])) create_args[[name]] <- optional[[name]]
}

t_create <- Sys.time()
rctd <- run_one(do.call(spacexr::create.RCTD, create_args))
if (inherits(rctd, "error")) {
  fail_rctd("spacexr::create.RCTD 失败", conditionMessage(rctd))
}
t_created <- Sys.time()

rctd <- run_one(spacexr::run.RCTD(rctd, doublet_mode = mode))
if (inherits(rctd, "error")) {
  fail_rctd(paste0("spacexr::run.RCTD(doublet_mode=", mode, ") 失败"), conditionMessage(rctd))
}
t_fit <- Sys.time()

# ---------------------------------------------------------------------------
# 导出（weights / results_df）
# ---------------------------------------------------------------------------

results <- rctd@results
used_spots <- colnames(rctd@spatialRNA@counts)
dropped <- setdiff(colnames(spot_mat), used_spots)
if (length(dropped)) {
  add_warn(sprintf("create.RCTD 按 UMI_min/counts_MIN 过滤掉 %d 个 spot（结果里没有这些 spot）",
                   length(dropped)))
}

raw_weights <- NULL
results_df <- NULL
if (!is.null(results$weights)) {
  raw_weights <- as.matrix(results$weights)
  results_df <- results$results_df
} else if (is.list(results) && length(results)) {
  # multi 模式：results 是每个 spot 一个 list（all_weights 为全部类型的权重）
  type_names <- rctd@cell_type_info$renorm[[2]]
  raw_weights <- matrix(0, nrow = length(results), ncol = length(type_names),
                        dimnames = list(used_spots, type_names))
  n_types_vec <- integer(length(results))
  min_score_vec <- rep(NA_real_, length(results))
  type_lists <- character(length(results))
  for (i in seq_along(results)) {
    entry <- results[[i]]
    vals <- entry$all_weights
    if (!is.null(vals)) {
      vals <- as.numeric(vals[type_names])
      vals[!is.finite(vals)] <- 0
      raw_weights[i, ] <- vals
    }
    lst <- entry$cell_type_list
    n_types_vec[i] <- if (is.null(lst)) 0L else length(lst)
    type_lists[i] <- if (is.null(lst)) "" else paste(lst, collapse = "+")
    if (!is.null(entry$min_score)) min_score_vec[i] <- as.numeric(entry$min_score)
  }
  results_df <- data.frame(
    spot_class = rep("multi", length(results)),
    n_types = n_types_vec,
    types = type_lists,
    min_score = min_score_vec,
    row.names = used_spots,
    stringsAsFactors = FALSE
  )
  add_warn("multi 模式：RCTD 的 results 是逐 spot 列表，weights 由各 spot 的 all_weights 汇总（types 列记录入选类型组合）")
} else {
  fail_runtime(sprintf("spacexr run.RCTD 的 results 结构无法识别（mode=%s）：%s",
                       mode, paste(names(results), collapse = ",")))
}

if (is.null(raw_weights) || nrow(raw_weights) == 0L) {
  fail_runtime("RCTD 没有返回任何 spot 的权重（可能所有 spot 都被 UMI 过滤掉了）")
}

cell_type_names <- colnames(raw_weights)
if (is.null(cell_type_names) || any(!nzchar(cell_type_names))) {
  fail_runtime("RCTD 权重矩阵没有细胞类型列名")
}

# 行归一化：spacexr 官方 normalize_weights()（full 模式权重只保证 ≤1）
row_sums <- rowSums(raw_weights)
n_zero <- sum(row_sums <= 0)
weights_norm <- raw_weights
if (any(row_sums > 0)) {
  weights_norm[row_sums > 0, ] <- spacexr::normalize_weights(raw_weights[row_sums > 0, , drop = FALSE])
}
if (n_zero > 0) {
  weights_norm[row_sums <= 0, ] <- 0
  add_warn(sprintf("%d 个 spot 的 RCTD 权重全为 0（无法解释的 spot）：比例记为 0", n_zero))
}

write.csv(as.matrix(weights_norm), file.path(out_dir, "tables", "weights.csv"),
          quote = TRUE)
write.csv(as.matrix(raw_weights), file.path(out_dir, "tables", "weights_raw.csv"),
          quote = TRUE)
if (!is.null(results_df)) {
  write.csv(results_df, file.path(out_dir, "tables", "results_df.csv"),
            quote = TRUE, row.names = TRUE)
} else {
  # full 模式没有 spot 分类：写空表（表头在），不伪造 spot_class
  write.csv(data.frame(spot_class = character(0), stringsAsFactors = FALSE),
            file.path(out_dir, "tables", "results_df.csv"), quote = TRUE,
            row.names = FALSE)
  add_warn("full 模式不产生 spot 分类（spacexr 的 results_df 只在 doublet/multi 模式生成）：results_df.csv 为空表")
}

n_dropped <- length(dropped)
result <- list(
  method = "RCTD",
  mode = mode,
  cell_types = cell_type_names,
  n_spots = length(used_spots),
  n_cells = ncol(ref_mat),
  n_genes = length(common),
  n_spots_input = ncol(spot_mat),
  n_spots_dropped = n_dropped,
  n_cells_per_type = as.list(as.integer(type_counts)),
  cell_type_col = cell_type_col,
  coords = coords_kind,
  max_cores = max_cores,
  cell_min = cell_min,
  weights_normalized = TRUE,
  raw_rowsum_median = if (length(row_sums)) as.numeric(median(row_sums)) else NA_real_,
  raw_rowsum_min = if (length(row_sums)) as.numeric(min(row_sums)) else NA_real_,
  raw_rowsum_max = if (length(row_sums)) as.numeric(max(row_sums)) else NA_real_,
  has_results_df = !is.null(results_df),
  spacexr_version = as.character(utils::packageVersion("spacexr")),
  r_version = R.version.string,
  elapsed_s = as.numeric(difftime(t_fit, t_start, units = "secs")),
  fit_s = as.numeric(difftime(t_fit, t_created, units = "secs")),
  warnings = as.list(WARNINGS)
)
names(result$n_cells_per_type) <- cell_types

tryCatch(
  write_json(result, file.path(out_dir, "result.json")),
  error = function(e) fail_runtime(paste0("result.json 写出失败：", conditionMessage(e)))
)

cat(sprintf("RCTD done: mode=%s, %d spots × %d types, %.1fs\n",
            mode, length(used_spots), length(cell_type_names), result$elapsed_s),
    file = stderr())
quit(save = "no", status = 0)
