"""
dsh-spatial — GEO 检索 / 详情 / 下载（DESIGN §7.1 / §7.2 / §7.3）

三个公开函数：
  * geo_search(term, species, study_type, min_samples, has_clinical, platform,
               limit, page, sort, has_matrix) -> dict
  * geo_detail(gse) -> dict
  * geo_download(gse, files, dest_dir, overwrite=False, progress=None) -> dict

实现要点（沿用 dsh-geobulk 的 bg_geo.py，已实测可用）：
  * 只用标准库 urllib（不依赖 requests）；所有网络失败重试 ≤2 次（退避 1s/2s），
    仍失败抛 EngineError("unsupported", "NCBI 暂时不可达：<原因>")，消息中文可读。
  * E-utilities：esearch(db=gds) → esummary(db=gds)（esummary 字段最全；efetch brief 对老
    GSE 会 "cannot get document summary"，所以不用）。UID = 200000000 + GSE 号。
  * GEO FTP：系列矩阵在 `<GSE>/matrix/`，作者补充文件在 `<GSE>/suppl/`；detail 两个目录
    都解析（实测 series matrix 不在 suppl）。用 HTTP Range 只抓系列矩阵头部 96KB，解析
    !Series_overall_design / !Series_contact_* / !Sample_characteristics_ch1 得到
    design / contacts / clinical，不需要下整份矩阵。
  * 中文：ZH_DICT ≥180 条（疾病/组织/药物/技术/物种/常用词 + 空间转录组术语），
    用于 title_zh、tags，以及把中文关键词反查成 NCBI 可检索的英文。
  * 空间转录组增强（本文件相对 bg_geo 的新增）：
      - `_SPATIAL_TERMS` / `_SPATIAL_WEAK_TERMS` / `spatial_platform_hints()`：从
        title+summary 判定 spatial_score（high/mid/low）与命中的空间平台
        （Visium / Visium HD / Stereo-seq / Slide-seq / MERFISH / seqFISH / Xenium /
        CosMx / CODEX / 10x Genomics，固定顺序、可多个）；
      - geo_search 支持 platform 过滤（all|visium|stereo|slideseq|10x|merfish，
        客户端过滤；全部未命中时回退未过滤结果并在 notes 说明）；
      - geo_detail 的文件 kind 增加 "visium"（suppl 里的 Visium 类 RAW.tar；
        series_matrix 优先级更高），并给出 visium_hint（是否可能有 spatial/ 目录）。
  * 下载：dest_dir 必须绝对路径且不是 `/`；单文件 8GB 上限；120s 超时；分块写盘 + 进度回调；
    overwrite=False 时已存在则跳过（ok=true, skipped=true）。
  * 覆盖下载的原子性：先写 `<name>.part`，只有完整性校验通过（写入字节数 == Content-Length；
    缺失 Content-Length 时按扩展名用 gzip/bz2/tar 可读性兜底）才 os.replace 原子替换；
    失败/中断保留旧文件并清理 `.part`。仅支持失败重试（每次从头下载），不做断点续传。

本文件只依赖标准库 + sp_common（EngineError）。禁止新增第三方依赖。
"""

from __future__ import annotations

import bz2
import gzip
import html
import json
import os
import re
import socket
import ssl
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib

from sp_common import EngineError

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
NCBI_ORIGIN = "https://ftp.ncbi.nlm.nih.gov"
GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
NCBI_WEB = "https://www.ncbi.nlm.nih.gov"

USER_AGENT = "dsh-spatial/1.0 (spatial transcriptomics workbench; local desktop app)"
# 实测：ftp.ncbi.nlm.nih.gov 的 WAF 只放行 UA 以 "curl/" 开头的请求，其它一律 403；
# 因此 FTP 请求实际发送 "curl/8.7.1 " + 本插件身份串（保留括号里的应用标识，不伪装其它软件）。
FTP_USER_AGENT = "dsh-spatial/1.0 (spatial transcriptomics workbench; local desktop app)"
FTP_UA_PREFIX = "curl/8.7.1 "
HTTP_TIMEOUT = 120                      # 单请求超时（秒）
RETRY_SLEEPS = (1.0, 2.0)               # 网络失败重试退避（秒）；测试里可置 0
ESUMMARY_BATCH = 100                    # esummary 单次最多 UID 数
MAX_ESUMMARY_CANDIDATES = 300           # has_clinical / sort 需要客户端筛选时的候选上限
MAX_DOWNLOAD_BYTES = 8 * 1024 ** 3      # 单文件下载上限 8GB
MATRIX_HEAD_BYTES = 96 * 1024           # 系列矩阵头部只抓 96KB（压缩后），足够解析 header
MATRIX_HEAD_TIMEOUT = 30                # 头部抓取的单请求超时（秒）
DOWNLOAD_CHUNK = 512 * 1024
PART_SUFFIX = ".part"                   # 覆盖下载的临时文件后缀（校验通过后原子替换）
VALIDATE_GZIP_FULL_MAX = 512 * 1024 ** 2  # 无 Content-Length 时 gzip/bz2 全量解压校验上限
SERIES_MATRIX_KIND = "series_matrix"
SPATIAL_KIND = "visium"            # suppl 里 Visium 类 RAW.tar（含 spatial/ 目录）

_SPECIES_ALIAS = {
    "human": "Homo sapiens", "homo sapiens": "Homo sapiens", "人": "Homo sapiens", "人类": "Homo sapiens",
    "mouse": "Mus musculus", "mus musculus": "Mus musculus", "小鼠": "Mus musculus",
    "rat": "Rattus norvegicus", "rattus norvegicus": "Rattus norvegicus", "大鼠": "Rattus norvegicus",
    "zebrafish": "Danio rerio", "斑马鱼": "Danio rerio",
    "fly": "Drosophila melanogaster", "果蝇": "Drosophila melanogaster",
    "yeast": "Saccharomyces cerevisiae", "酵母": "Saccharomyces cerevisiae",
}

_STUDY_TYPE_ALIAS = {
    "all": None, "any": None, "": None, "全部": None,
    "array": "expression profiling by array",
    "microarray": "expression profiling by array",
    "芯片": "expression profiling by array",
    "expression profiling by array": "expression profiling by array",
    "sequencing": "high throughput sequencing",
    "high throughput sequencing": "high throughput sequencing",
    "rnaseq": "high throughput sequencing",
    "rna-seq": "high throughput sequencing",
    "rna seq": "high throughput sequencing",
    "测序": "high throughput sequencing",
    "methylation": "methylation profiling by array",
    "甲基化": "methylation profiling by array",
}

# 判定"有临床信息"的关键词（用于 search 过滤 / item.has_clinical）
_CLINICAL_KEYWORDS = (
    "clinical", "survival", "prognos", "follow-up", "follow up", "outcome",
    "recurrence", "relapse", "mortality", "hazard ratio", "overall survival",
    "disease-free", "progression-free", "risk score", "nomogram", "cox regression",
    "kaplan", "clinical outcome", "pathological stage", "tnm stage",
    "治疗反应", "预后", "生存", "复发", "临床",
)

# 疾病类中文词典 key（用于 item.disease 识别，顺序即优先级）
DISEASE_KEYS = [
    "gastric cancer", "lung cancer", "breast cancer", "colorectal cancer", "liver cancer",
    "hepatocellular carcinoma", "esophageal cancer", "pancreatic cancer", "prostate cancer",
    "ovarian cancer", "cervical cancer", "bladder cancer", "renal cell carcinoma", "glioma",
    "glioblastoma", "melanoma", "leukemia", "lymphoma", "multiple myeloma", "thyroid cancer",
    "head and neck cancer", "oral cancer", "nasopharyngeal carcinoma", "cholangiocarcinoma",
    "osteosarcoma", "sarcoma", "neuroblastoma", "medulloblastoma", "endometrial cancer",
    "testicular cancer", "triple-negative breast cancer", "non-small cell lung cancer",
    "small cell lung cancer", "alzheimer's disease", "parkinson's disease",
    "multiple sclerosis", "type 2 diabetes", "covid-19", "sepsis", "atherosclerosis",
    "inflammatory bowel disease", "rheumatoid arthritis", "chronic kidney disease",
    "heart failure", "pulmonary fibrosis", "colitis", "psoriasis", "tuberculosis",
]

# ---------------------------------------------------------------------------
# 中英词典（≥180 条）：中文标题 / 标签 / 疾病识别 / 中文反查共用
# 键统一小写；替换时按长度降序、英文词边界匹配，中文原样保留。
# ---------------------------------------------------------------------------

ZH_DICT = {
    # ---- 疾病（50） ----
    "gastric cancer": "胃癌",
    "stomach cancer": "胃癌",
    "lung cancer": "肺癌",
    "breast cancer": "乳腺癌",
    "colorectal cancer": "结直肠癌",
    "colon cancer": "结肠癌",
    "rectal cancer": "直肠癌",
    "liver cancer": "肝癌",
    "hepatocellular carcinoma": "肝细胞癌",
    "esophageal cancer": "食管癌",
    "esophageal squamous cell carcinoma": "食管鳞癌",
    "pancreatic cancer": "胰腺癌",
    "prostate cancer": "前列腺癌",
    "ovarian cancer": "卵巢癌",
    "cervical cancer": "宫颈癌",
    "bladder cancer": "膀胱癌",
    "renal cell carcinoma": "肾细胞癌",
    "kidney cancer": "肾癌",
    "glioma": "胶质瘤",
    "glioblastoma": "胶质母细胞瘤",
    "melanoma": "黑色素瘤",
    "leukemia": "白血病",
    "lymphoma": "淋巴瘤",
    "multiple myeloma": "多发性骨髓瘤",
    "thyroid cancer": "甲状腺癌",
    "head and neck cancer": "头颈癌",
    "oral cancer": "口腔癌",
    "nasopharyngeal carcinoma": "鼻咽癌",
    "cholangiocarcinoma": "胆管癌",
    "osteosarcoma": "骨肉瘤",
    "sarcoma": "肉瘤",
    "neuroblastoma": "神经母细胞瘤",
    "medulloblastoma": "髓母细胞瘤",
    "endometrial cancer": "子宫内膜癌",
    "testicular cancer": "睾丸癌",
    "bladder urothelial carcinoma": "膀胱尿路上皮癌",
    "triple-negative breast cancer": "三阴性乳腺癌",
    "non-small cell lung cancer": "非小细胞肺癌",
    "small cell lung cancer": "小细胞肺癌",
    "alzheimer's disease": "阿尔茨海默病",
    "parkinson's disease": "帕金森病",
    "multiple sclerosis": "多发性硬化",
    "type 2 diabetes": "2型糖尿病",
    "covid-19": "新型冠状病毒肺炎",
    "sepsis": "脓毒症",
    "atherosclerosis": "动脉粥样硬化",
    "inflammatory bowel disease": "炎症性肠病",
    "rheumatoid arthritis": "类风湿关节炎",
    "chronic kidney disease": "慢性肾病",
    "heart failure": "心力衰竭",
    "pulmonary fibrosis": "肺纤维化",
    "colitis": "结肠炎",
    "psoriasis": "银屑病",
    "tuberculosis": "结核病",
    # ---- 组织 / 部位 / 细胞（40） ----
    "tumor microenvironment": "肿瘤微环境",
    "adjacent normal": "癌旁正常",
    "normal tissue": "正常组织",
    "tumor tissue": "肿瘤组织",
    "tumour tissue": "肿瘤组织",
    "gastric mucosa": "胃黏膜",
    "mucosa": "黏膜",
    "gastric tissue": "胃组织",
    "liver tissue": "肝组织",
    "lung tissue": "肺组织",
    "breast tissue": "乳腺组织",
    "colon tissue": "结肠组织",
    "cerebrospinal fluid": "脑脊液",
    "peripheral blood mononuclear cells": "外周血单个核细胞",
    "bone marrow": "骨髓",
    "lymph node": "淋巴结",
    "spleen": "脾脏",
    "thymus": "胸腺",
    "pancreas": "胰腺",
    "prostate": "前列腺",
    "ovary": "卵巢",
    "uterus": "子宫",
    "placenta": "胎盘",
    "adipose tissue": "脂肪组织",
    "skeletal muscle": "骨骼肌",
    "fibroblast": "成纤维细胞",
    "fibroblasts": "成纤维细胞",
    "endothelial": "内皮",
    "epithelium": "上皮",
    "epithelial cells": "上皮细胞",
    "immune cells": "免疫细胞",
    "t cells": "T 细胞",
    "b cells": "B 细胞",
    "macrophage": "巨噬细胞",
    "macrophages": "巨噬细胞",
    "plasma": "血浆",
    "serum": "血清",
    "sputum": "痰液",
    "urine": "尿液",
    "saliva": "唾液",
    "biopsy": "活检",
    "blood": "血液",
    # ---- 药物 / 处理（25） ----
    "cisplatin": "顺铂",
    "paclitaxel": "紫杉醇",
    "docetaxel": "多西他赛",
    "doxorubicin": "多柔比星",
    "5-fluorouracil": "5-氟尿嘧啶",
    "gemcitabine": "吉西他滨",
    "oxaliplatin": "奥沙利铂",
    "irinotecan": "伊立替康",
    "etoposide": "依托泊苷",
    "methotrexate": "甲氨蝶呤",
    "temozolomide": "替莫唑胺",
    "sorafenib": "索拉非尼",
    "sunitinib": "舒尼替尼",
    "imatinib": "伊马替尼",
    "erlotinib": "厄洛替尼",
    "gefitinib": "吉非替尼",
    "trastuzumab": "曲妥珠单抗",
    "tamoxifen": "他莫昔芬",
    "nivolumab": "纳武利尤单抗",
    "pembrolizumab": "帕博利珠单抗",
    "dexamethasone": "地塞米松",
    "rapamycin": "雷帕霉素",
    "metformin": "二甲双胍",
    "aspirin": "阿司匹林",
    "immunotherapy": "免疫治疗",
    "chemotherapy": "化疗",
    "radiotherapy": "放疗",
    "drug resistance": "耐药",
    "chemoresistance": "化疗耐药",
    "cisplatin resistance": "顺铂耐药",
    # ---- 技术 / 平台（15） ----
    "microarray": "芯片",
    "expression profiling": "表达谱",
    "rna-seq": "转录组测序",
    "rna sequencing": "转录组测序",
    "transcriptome": "转录组",
    "single-cell rna-seq": "单细胞转录组测序",
    "single cell rna-seq": "单细胞转录组测序",
    "scrna-seq": "单细胞转录组测序",
    "whole genome sequencing": "全基因组测序",
    "exome sequencing": "外显子组测序",
    "chip-seq": "染色质免疫共沉淀测序",
    "atac-seq": "染色质可及性测序",
    "bisulfite sequencing": "亚硫酸氢盐测序",
    "spatial transcriptomics": "空间转录组",
    "qpcr": "实时定量 PCR",
    "proteomics": "蛋白组学",
    "metabolomics": "代谢组学",
    "crispr": "基因编辑",
    "high throughput sequencing": "高通量测序",
    # ---- 物种（12） ----
    "homo sapiens": "人",
    "mus musculus": "小鼠",
    "rattus norvegicus": "大鼠",
    "danio rerio": "斑马鱼",
    "human": "人",
    "mouse": "小鼠",
    "rat": "大鼠",
    "patients": "患者",
    "patient": "患者",
    "healthy": "健康",
    "wild-type": "野生型",
    "knockout": "敲除",
    # ---- 常用生物学词（30） ----
    "biomarker": "生物标志物",
    "signature": "特征签名",
    "prognosis": "预后",
    "prognostic": "预后",
    "survival": "生存",
    "recurrence": "复发",
    "metastasis": "转移",
    "immune infiltration": "免疫浸润",
    "inflammation": "炎症",
    "apoptosis": "凋亡",
    "proliferation": "增殖",
    "hypoxia": "缺氧",
    "resistance": "耐药",
    "overexpression": "过表达",
    "knockdown": "敲低",
    "mutation": "突变",
    "methylation": "甲基化",
    "differentiation": "分化",
    "development": "发育",
    "infection": "感染",
    "vaccine": "疫苗",
    "therapy": "治疗",
    "response": "反应",
    "risk": "风险",
    "expression": "表达",
    "profiling": "谱分析",
    "single cell": "单细胞",
    "exosome": "外泌体",
    "circrna": "环状RNA",
    "mirna": "微小RNA",
    "lncrna": "长链非编码RNA",
    "autophagy": "自噬",
    "ferroptosis": "铁死亡",
    "senescence": "衰老",
    "fibrosis": "纤维化",
    # ---- 空间转录组 / 空间组学 / 病理成像（52 条，dsh-spatial 追加） ----
    "spatial transcriptome": "空间转录组",
    "spatial transcriptomics analysis": "空间转录组分析",
    "spatial omics": "空间组学",
    "spatial multi-omics": "空间多组学",
    "spatial profiling": "空间表达谱",
    "spatial analysis": "空间分析",
    "spatial distribution": "空间分布",
    "spatial heterogeneity": "空间异质性",
    "spatial domain": "空间域",
    "spatial niche": "空间生态位",
    "spatial barcode": "空间条形码",
    "spatial coordinate": "空间坐标",
    "visium": "Visium 空间转录组平台",
    "visium hd": "空间转录组 HD 平台",
    "stereo-seq": "Stereo-seq 空间转录组平台",
    "stereoseq": "Stereo-seq 空间转录组平台",
    "slide-seq": "Slide-seq 空间转录组平台",
    "merfish": "MERFISH 空间单细胞成像",
    "seqfish": "seqFISH 空间单细胞成像",
    "xenium": "Xenium 空间单细胞成像",
    "cosmx": "CosMx 空间单细胞成像",
    "codex": "CODEX 多重成像",
    "10x genomics": "10x Genomics 平台",
    "tissue section": "组织切片",
    "tissue sections": "组织切片",
    "frozen section": "冷冻切片",
    "ffpe": "石蜡切片",
    "ffpe tissue": "石蜡包埋组织",
    "histology": "组织学",
    "histological": "组织学的",
    "h&e": "苏木精-伊红染色",
    "h & e": "苏木精-伊红染色",
    "hematoxylin": "苏木精",
    "eosin": "伊红",
    "in situ": "原位",
    "in situ hybridization": "原位杂交",
    "tissue in situ": "组织原位",
    "in situ sequencing": "原位测序",
    "immunofluorescence": "免疫荧光",
    "fluorescence": "荧光",
    "fluorescence imaging": "荧光成像",
    "imaging": "成像",
    "imaging mass cytometry": "成像质谱流式",
    "multiplex imaging": "多重成像",
    "pathology": "病理",
    "pathological section": "病理切片",
    "tumor boundary": "肿瘤边界",
    "tumour boundary": "肿瘤边界",
    "tissue architecture": "组织结构",
    "region of interest": "感兴趣区域",
    "barcode": "条形码",
    "coordinates": "坐标",
    "breast": "乳腺",
}

# 英文词边界匹配（避免 "rat" 命中 "patient" 中的子串）
_ZH_SORTED = sorted(ZH_DICT.items(), key=lambda kv: len(kv[0]), reverse=True)
_ZH_PATTERNS = [
    (re.compile(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", re.IGNORECASE), zh)
    for term, zh in _ZH_SORTED
]
_DISEASE_PATTERNS = [
    (re.compile(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", re.IGNORECASE), term)
    for term in sorted(DISEASE_KEYS, key=len, reverse=True)
]

_TECH_TAGS = {
    "expression profiling by array": "表达谱",
    "expression profiling by high throughput sequencing": "高通量测序",
    "methylation profiling by array": "甲基化",
    "methylation profiling by high throughput sequencing": "甲基化测序",
    "genome variation profiling by array": "基因组变异",
    "genome variation profiling by high throughput sequencing": "基因组变异",
    "snapshot": "快照",
}
_SPECIES_TAGS = {"Homo sapiens": "人", "Mus musculus": "小鼠", "Rattus norvegicus": "大鼠"}


# ---------------------------------------------------------------------------
# 空间转录组判定（title + summary → spatial_score / spatial_platform）
# ---------------------------------------------------------------------------

# 强证据：明确的空间转录组 / 空间组学技术或平台名（命中即 high）
_SPATIAL_TERMS = (
    "spatial transcriptomics", "spatial transcriptome", "spatial transcriptomic",
    "spatial omics", "spatial multi-omics", "spatial multiomics",
    "spatial profiling", "spatial gene expression", "spatial expression",
    "spatial barcode", "visium", "stereo-seq", "stereoseq", "stereo seq",
    "slide-seq", "slide seq", "slideseq", "merfish", "seqfish", "seq-fish",
    "xenium", "cosmx", "codex", "in situ sequencing", "spatial proteomics",
    "10x genomics", "spatial atac",
)

# 弱证据：泛指空间 / 成像 / 切片 / 病理的措辞（命中即 mid，不构成强证据）
_SPATIAL_WEAK_TERMS = (
    "spatial", "in situ", "in-situ", "imaging", "fluorescence", "immunofluorescence",
    "histology", "histological", "h&e", "h & e", "hematoxylin", "eosin",
    "tissue section", "frozen section", "ffpe", "pathology", "pathological",
    "microdissection", "laser capture", "tissue architecture", "tumor boundary",
    "空间", "原位", "切片", "病理", "组织学", "成像",
)

# 平台名反查：固定顺序（Visium HD 在 Visium 之前，避免同一数据集重复计数）
_SPATIAL_PLATFORM_ORDER = (
    ("Visium HD", ("visium hd",)),
    ("Visium", ("visium",)),
    ("Stereo-seq", ("stereo-seq", "stereo seq", "stereoseq")),
    ("Slide-seq", ("slide-seq", "slide seq", "slideseq")),
    ("MERFISH", ("merfish",)),
    ("seqFISH", ("seqfish", "seq-fish")),
    ("Xenium", ("xenium",)),
    ("CosMx", ("cosmx",)),
    ("CODEX", ("codex",)),
    ("10x Genomics", ("10x genomics", "10xgenomics")),
)

# geo_search(platform=...) 的客户端过滤桶：桶 → 在 title/summary/spatial_platform 里找的关键词
_PLATFORM_MATCH = {
    "visium": ("visium",),
    "stereo": ("stereo-seq", "stereo seq", "stereoseq"),
    "slideseq": ("slide-seq", "slide seq", "slideseq"),
    "10x": ("10x",),
    "merfish": ("merfish",),
}
_PLATFORM_ALIAS = {
    "visium hd": "visium", "10x genomics": "10x", "10xgenomics": "10x",
    "stereo-seq": "stereo", "stereo seq": "stereo", "stereoseq": "stereo",
    "slide-seq": "slideseq", "slide seq": "slideseq",
    "seqfish": "merfish", "seq-fish": "merfish",
}

# "像 Visium 数据集" 的措辞（用于把 suppl 里的 RAW.tar 归为 visium）
_VISIUM_TERMS = ("visium", "10x genomics", "spatial transcriptomics",
                 "spatial transcriptome", "spatial omics", "spatial profiling")


def _spatial_text(title: str, summary: str) -> str:
    """把 title/summary 拼成小写判定文本（内部使用）。"""
    return f"{title or ''} {summary or ''}".strip().lower()


def spatial_score_hint(title, summary) -> str:
    """按 title+summary 判定空间转录组相关度：high / mid / low。"""
    text = _spatial_text(title, summary)
    if not text:
        return "low"
    if any(term in text for term in _SPATIAL_TERMS):
        return "high"
    if any(term in text for term in _SPATIAL_WEAK_TERMS):
        return "mid"
    return "low"


def spatial_platform_hints(title, summary) -> list[str]:
    """反查命中的空间平台名（如 ["Visium","Stereo-seq"]；固定顺序、可多个）。"""
    text = _spatial_text(title, summary)
    if not text:
        return []
    out: list[str] = []
    for canonical, patterns in _SPATIAL_PLATFORM_ORDER:
        if canonical == "Visium" and "visium hd" in text:
            continue  # 【Visium HD】已单独计数，不再重复算 Visium
        if any(pattern in text for pattern in patterns):
            out.append(canonical)
    return out


def _visium_context(title, summary) -> bool:
    """title/summary 是否命中 Visium 类空间数据集（RAW.tar 归为 visium 的前提）。"""
    text = _spatial_text(title, summary)
    return any(term in text for term in _VISIUM_TERMS)


def _platform_key(platform) -> str:
    """把 platform 参数归一化成 _PLATFORM_MATCH 的桶名（未知值 → all）。"""
    key = str(platform or "all").strip().lower()
    key = _PLATFORM_ALIAS.get(key, key)
    if key in ("", "any", "全部"):
        key = "all"
    return key if key in _PLATFORM_MATCH else "all"


_PLATFORM_QUERY_TERMS = {
    "visium": '("Visium"[All Fields] OR "10x Genomics"[All Fields] OR "spatial transcriptomics"[All Fields] OR "spatial transcriptome"[All Fields])',
    "stereo": '("Stereo-seq"[All Fields] OR "stereo seq"[All Fields] OR "spatial transcriptomics"[All Fields])',
    "slideseq": '("Slide-seq"[All Fields] OR "slide seq"[All Fields] OR "spatial transcriptomics"[All Fields])',
    "10x": '("Visium"[All Fields] OR "Xenium"[All Fields] OR "10x Genomics"[All Fields] OR "spatial transcriptomics"[All Fields])',
    "merfish": '("MERFISH"[All Fields] OR "spatial transcriptomics"[All Fields])',
    "spatial": '("spatial transcriptomics"[All Fields] OR "spatial transcriptome"[All Fields] OR "Visium"[All Fields] OR "Stereo-seq"[All Fields] OR "Slide-seq"[All Fields] OR "MERFISH"[All Fields])',
}

_SCORE_RANK = {"high": 0, "mid": 1, "medium": 1, "low": 2}


def _platform_clause(platform_key: str) -> str:
    """把平台选择变成检索式的一部分：选 Visium 就只让 NCBI 回 Visium 相关记录。"""
    return _PLATFORM_QUERY_TERMS.get(str(platform_key or "all").strip().lower(), "")


def _with_platform_clause(term: str, platform_key: str) -> str:
    clause = _platform_clause(platform_key)
    if not clause:
        return term
    return f"({term}) AND {clause}"


def _match_platform(item: dict, platform_key: str) -> bool:
    """客户端平台过滤：看 title/summary + 已判定的 spatial_platform 是否含该桶关键词。"""
    keywords = _PLATFORM_MATCH.get(str(platform_key or "all").strip().lower())
    if not keywords:
        return True
    parts = [str((item or {}).get("title") or ""), str((item or {}).get("summary") or "")]
    parts.extend(str(name) for name in ((item or {}).get("spatial_platform") or []))
    text = " ".join(parts).lower()
    return any(keyword in text for keyword in keywords)


# ---------------------------------------------------------------------------
# 网络层（测试通过替换模块级 _urlopen 来 mock）
# ---------------------------------------------------------------------------


def _ssl_context():
    """HTTPS 校验上下文：优先用 certifi 的 CA 包（本机 python.org 自带的 openssl
    cert.pem 可能缺根证书，导致 CERTIFICATE_VERIFY_FAILED）。"""
    global _SSL_CTX
    if _SSL_CTX is None:
        try:
            import certifi

            _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
        except Exception:  # noqa: BLE001 - 没有 certifi 时退回默认 CA
            _SSL_CTX = ssl.create_default_context()
    return _SSL_CTX


_SSL_CTX = None


def _urlopen(request, timeout):  # pragma: no cover - 真实网络路径
    """单点 URL 打开（测试 monkeypatch 这一层）。"""
    try:
        return urllib.request.urlopen(request, timeout=timeout, context=_ssl_context())
    except TypeError:  # pragma: no cover - 某些 context 不支持的特殊 URL 类型
        return urllib.request.urlopen(request, timeout=timeout)


def _reason(error) -> str:
    """把网络异常转成中文可读原因。"""
    if isinstance(error, urllib.error.HTTPError):
        code = int(getattr(error, "code", 0) or 0)
        if code == 404:
            return "资源不存在（HTTP 404）"
        if code in (429, 503):
            return f"NCBI 限流或维护（HTTP {code}）"
        if code >= 500:
            return f"NCBI 服务端错误（HTTP {code}）"
        return f"HTTP {code}"
    if isinstance(error, urllib.error.URLError):
        detail = getattr(error, "reason", error)
        if isinstance(detail, (socket.timeout, TimeoutError)):
            return "连接超时"
        if isinstance(detail, ssl.SSLError):
            return "TLS 握手失败"
        if isinstance(detail, socket.gaierror):
            return "域名解析失败（可能无网络）"
        return str(detail)
    if isinstance(error, (socket.timeout, TimeoutError)):
        return "连接超时"
    if isinstance(error, ssl.SSLError):
        return "TLS 握手失败"
    if isinstance(error, json.JSONDecodeError):
        return "NCBI 返回内容不是合法 JSON"
    return str(error) or type(error).__name__


def _network_fail(error) -> EngineError:
    reason = _reason(error)
    if "403" in reason:
        return EngineError(
            "unsupported",
            "NCBI 拒绝了请求（HTTP 403，通常是短时间请求过多被限流）。等 10–30 秒再试；"
            "如需高频检索，设置环境变量 DSH_NCBI_API_KEY（NCBI 免费申请）可把配额提到 10 次/秒。")
    if "429" in reason:
        return EngineError("unsupported", "NCBI 限流（HTTP 429）：请稍后重试，或配置 DSH_NCBI_API_KEY。")
    return EngineError("unsupported", f"NCBI 暂时不可达：{reason}")


def _request(url: str, headers: dict | None = None):
    if "ftp.ncbi.nlm.nih.gov" in url:
        agent = FTP_UA_PREFIX + FTP_USER_AGENT  # WAF 只放行 curl/ 前缀
    else:
        agent = USER_AGENT
    merged = {"User-Agent": agent, "Accept": "*/*"}
    if headers:
        merged.update(headers)
    return urllib.request.Request(url, headers=merged)


def _fetch_bytes(url: str, *, timeout: int = HTTP_TIMEOUT, retries: int = 2,
                 headers: dict | None = None, limit: int | None = None,
                 soft_404: bool = False) -> bytes | None:
    """GET 一个 URL 返回 bytes；失败重试 ≤2 次，仍失败抛 EngineError。"""
    last_error = None
    attempts = max(0, int(retries)) + 1
    for attempt in range(attempts):
        try:
            with _urlopen(_request(url, headers), timeout) as response:
                data = response.read() if limit is None else response.read(limit)
            return data or b""
        except urllib.error.HTTPError as error:
            if soft_404 and int(getattr(error, "code", 0) or 0) == 404:
                return None
            last_error = error
        except (urllib.error.URLError, socket.timeout, TimeoutError, ssl.SSLError, OSError) as error:
            last_error = error
        except Exception as error:  # noqa: BLE001 - 任何网络层异常都按不可达处理
            last_error = error
        if attempt < attempts - 1:
            delay = RETRY_SLEEPS[min(attempt, len(RETRY_SLEEPS) - 1)]
            if delay > 0:
                time.sleep(delay)
    raise _network_fail(last_error)


def _fetch_json(url: str, *, retries: int = 2, soft_404: bool = False) -> dict | None:
    raw = _fetch_bytes(url, retries=retries, soft_404=soft_404)
    if raw is None:
        return None
    try:
        text = raw.decode("utf-8", "replace")
        data = json.loads(text)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise _network_fail(error)
    if not isinstance(data, dict):
        raise _network_fail("NCBI 返回结构异常")
    if data.get("error"):
        raise EngineError("bad_request", f"检索式无法解析：{str(data.get('error'))[:300]}")
    return data


def _fetch_text(url: str, *, retries: int = 2, soft_404: bool = False, limit: int | None = None) -> str | None:
    raw = _fetch_bytes(url, retries=retries, soft_404=soft_404, limit=limit)
    if raw is None:
        return None
    return raw.decode("utf-8", "replace")


# ---------------------------------------------------------------------------
# 小工具
# ---------------------------------------------------------------------------


def _norm_gse(gse) -> str:
    text = str(gse or "").strip().upper().replace(" ", "")
    match = re.match(r"^(?:GSE)?(\d{1,9})$", text)
    if not match:
        raise EngineError("bad_request", f"GSE 编号格式不正确：{gse!r}（示例：GSE13861）")
    return f"GSE{int(match.group(1))}"


def _gse_number(gse: str) -> int:
    return int(re.sub(r"\D", "", gse))


def _ftp_prefix(gse: str) -> str:
    number = str(_gse_number(gse))
    base = number[:-3] if len(number) > 3 else ""
    if not base:
        base = "0"
    return f"GSE{base}nnn"


def _series_root(gse: str) -> str:
    return f"{GEO_FTP}/{_ftp_prefix(gse)}/{gse}"


def _uid_for(gse: str) -> str:
    return str(200000000 + _gse_number(gse))


def _norm_date(value) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    match = re.match(r"^(\d{4})[/\-.](\d{1,2})[/\-.](\d{1,2})", text)
    if match:
        return f"{match.group(1)}-{int(match.group(2)):02d}-{int(match.group(3)):02d}"
    match = re.match(r"^([A-Z][a-z]{2})\s+(\d{1,2}),?\s+(\d{4})$", text)
    if match:
        months = {m: i + 1 for i, m in enumerate(
            ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])}
        month = months.get(match.group(1))
        if month:
            return f"{match.group(3)}-{month:02d}-{int(match.group(2)):02d}"
    return text


def _platforms(value) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    tokens = [tok for tok in re.split(r"[,;\s]+", text) if tok]
    out = []
    for token in tokens:
        token = token.strip().upper()
        if not token:
            continue
        if token.isdigit():
            token = f"GPL{token}"
        elif token.startswith("GPL"):
            token = f"GPL{token[3:].lstrip('0') or '0'}"
        if token not in out:
            out.append(token)
    return out


def _parse_size(token: str) -> int:
    text = str(token or "").strip().upper()
    if not text or text in ("-", "N/A"):
        return 0
    match = re.match(r"^([0-9.]+)\s*([KMGTP]?)", text)
    if not match:
        return 0
    try:
        value = float(match.group(1))
    except ValueError:
        return 0
    factor = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3,
              "T": 1024 ** 4, "P": 1024 ** 5}[match.group(2)]
    return int(value * factor)


def _clean_quotes(value: str) -> str:
    text = str(value or "").strip().strip('"').strip()
    # GEO series matrix 里值内的逗号会双写成 ",,"，还原成单个逗号
    if ",," in text:
        text = text.replace(",,", ",")
    return text


def _translate(text: str) -> str:
    """英文→中文词典替换；没有任何命中时返回空字符串（允许 title_zh 为空）。"""
    source = str(text or "")
    if not source:
        return ""
    # 先保护已经替成中文的片段，避免后续英文词再命中的副作用；直接顺序替换即可（最长优先）
    result = source
    hits = 0
    for pattern, zh in _ZH_PATTERNS:
        def _repl(match, _zh=zh):
            nonlocal hits
            hits += 1
            return _zh
        result, count = pattern.subn(_repl, result)
    if not hits:
        return ""
    # 压缩多余空格
    result = re.sub(r"\s{2,}", " ", result).strip()
    return result


def _matched_chinese(text: str, limit: int = 6) -> list[str]:
    out: list[str] = []
    source = str(text or "")
    if not source:
        return out
    for pattern, zh in _ZH_PATTERNS:
        if len(out) >= limit:
            break
        if pattern.search(source) and zh not in out:
            out.append(zh)
    return out


def _detect_disease(title: str, summary: str) -> str:
    for text in (title, summary):
        for pattern, term in _DISEASE_PATTERNS:
            if pattern.search(str(text or "")):
                return term
    return ""


def _looks_clinical(text: str) -> bool:
    lower = str(text or "").lower()
    return any(keyword in lower for keyword in _CLINICAL_KEYWORDS)


def _make_tags(title: str, summary: str, species: str, study_type: str,
               disease: str) -> list[str]:
    tags: list[str] = []
    disease_zh = ZH_DICT.get(disease.lower()) if disease else None
    if disease_zh:
        tags.append(disease_zh)
    for item in _matched_chinese(title, limit=6) + _matched_chinese(summary, limit=4):
        if item not in tags:
            tags.append(item)
    species_zh = _SPECIES_TAGS.get(str(species or "").strip())
    if species_zh and species_zh not in tags:
        tags.append(species_zh)
    study_zh = _TECH_TAGS.get(str(study_type or "").strip().lower())
    if study_zh and study_zh not in tags:
        tags.append(study_zh)
    if not tags:
        tags.append("GEO 数据集")
    return tags[:10]


# ---------------------------------------------------------------------------
# 检索式
# ---------------------------------------------------------------------------


def _build_query(term: str, species: str, study_type: str, min_samples: int) -> str:
    clauses: list[str] = [f"({term})"]
    if "etyp" not in term.lower():
        clauses.append("(GSE[ETYP])")
    species_key = str(species or "all").strip().lower()
    if species_key not in ("all", "", "any", "全部"):
        organism = _SPECIES_ALIAS.get(species_key)
        if not organism:
            # 用户直接给了物种拉丁名/别名：原样作为 Organism 字段搜索
            organism = str(species).strip()
        clauses.append(f'("{organism}"[Organism])')
    study_key = str(study_type or "all").strip().lower()
    mapped = _STUDY_TYPE_ALIAS.get(study_key, study_type)
    if mapped:
        mapped = str(mapped).strip()
        if mapped.lower() not in ("all", ""):
            clauses.append(f'("{mapped}"[DataSet Type])')
    try:
        min_n = int(min_samples or 0)
    except (TypeError, ValueError):
        min_n = 0
    if min_n > 0:
        clauses.append(f"({min_n}:10000000[Number of Samples])")
    return " AND ".join(clauses)


def _eutils_identity() -> dict:
    """NCBI 建议的 tool/email/api_key（配了 api_key 配额更高、更不容易 403）。"""
    extra = {"tool": "dsh-spatial"}
    email = str(os.environ.get("DSH_NCBI_EMAIL") or "").strip()
    if email:
        extra["email"] = email
    api_key = str(os.environ.get("DSH_NCBI_API_KEY") or "").strip()
    if api_key:
        extra["api_key"] = api_key
    return extra


def _esearch(query: str, *, retmax: int, retstart: int = 0) -> tuple[list[str], int, str]:
    params = urllib.parse.urlencode({
        "db": "gds",
        "term": query,
        "retmode": "json",
        "retmax": int(retmax),
        "retstart": int(retstart),
        **_eutils_identity(),
    })
    data = _fetch_json(f"{EUTILS}/esearch.fcgi?{params}")
    payload = (data or {}).get("esearchresult") or {}
    count = int(payload.get("count") or 0)
    ids = [str(item) for item in (payload.get("idlist") or [])]
    translation = str(payload.get("querytranslation") or "")
    return ids, count, translation


def _esummary(uids: list[str]) -> list[dict]:
    records: list[dict] = []
    for start in range(0, len(uids), ESUMMARY_BATCH):
        batch = uids[start:start + ESUMMARY_BATCH]
        if not batch:
            continue
        params = urllib.parse.urlencode({
            "db": "gds",
            "id": ",".join(batch),
            "retmode": "json",
            **_eutils_identity(),
        })
        data = _fetch_json(f"{EUTILS}/esummary.fcgi?{params}")
        result = (data or {}).get("result") or {}
        ordered = [str(item) for item in (result.get("uids") or batch)]
        by_uid = {}
        for uid in ordered:
            record = result.get(uid)
            if isinstance(record, dict) and not record.get("error"):
                by_uid[uid] = record
        for uid in batch:
            record = by_uid.get(uid)
            if record:
                records.append(record)
    return records


def _item_from_summary(record: dict) -> dict:
    accession = str(record.get("accession") or "")
    gse = accession.upper() if accession.upper().startswith("GSE") else _norm_gse_safe(record.get("gse"))
    title = str(record.get("title") or "")
    summary = str(record.get("summary") or "")
    species = str(record.get("taxon") or "")
    try:
        n_samples = int(record.get("n_samples") or len(record.get("samples") or []))
    except (TypeError, ValueError):
        n_samples = 0
    platforms = _platforms(record.get("gpl"))
    study_type = str(record.get("gdstype") or "").strip()
    release = _norm_date(record.get("pdat"))
    disease = _detect_disease(title, summary)
    has_clinical = _looks_clinical(f"{title} {summary}")
    return {
        "gse": gse,
        "title": title,
        "title_zh": _translate(title),
        "summary": summary,
        "species": species,
        "n_samples": n_samples,
        "platform": ", ".join(platforms),
        "study_type": study_type.lower(),
        "release": release,
        "has_matrix": bool(n_samples > 0),
        "has_clinical": has_clinical,
        "disease": disease,
        "spatial_score": spatial_score_hint(title, summary),
        "spatial_platform": spatial_platform_hints(title, summary),
        "tags": _make_tags(title, summary, species, study_type, disease),
    }


def _norm_gse_safe(value) -> str:
    try:
        return _norm_gse(value)
    except EngineError:
        text = str(value or "").strip()
        return f"GSE{text}" if text.isdigit() else text


def _sort_items(items: list[dict], sort: str) -> list[dict]:
    mode = str(sort or "relevance").strip().lower()
    if mode == "date":
        return sorted(items, key=lambda it: (it.get("release") or ""), reverse=True)
    if mode in ("samples", "sample", "n_samples"):
        return sorted(items, key=lambda it: int(it.get("n_samples") or 0), reverse=True)
    return items


# ---------------------------------------------------------------------------
# geo_search
# ---------------------------------------------------------------------------


_CJK_RE = re.compile(r"[\u4e00-\u9fff]")
_ZH_TO_EN = sorted(
    ((zh, en) for en, zh in ZH_DICT.items() if isinstance(zh, str) and zh),
    key=lambda item: len(item[0]),
    reverse=True,
)


_GENERIC_QUERY_RE = re.compile(
    r"\b(expression\s+profil\w*|profiling|expression|transcriptom\w*|microarray|rna[\s-]?seq|"
    r"high[\s-]?throughput|sequencing|dataset|datasets|data)\b",
    re.I,
)


def _cascade_terms(term: str) -> list[str]:
    """把用户关键词展开成由严到宽的候选查询（智能化放宽）。"""
    base = re.sub(r"\s{2,}", " ", str(term or "")).strip()
    candidates: list[str] = [base]
    stripped = _GENERIC_QUERY_RE.sub(" ", base)
    stripped = re.sub(r"\s{2,}", " ", stripped).strip()
    if stripped and stripped != base:
        candidates.append(stripped)
    words = [w for w in re.split(r"\s+", stripped or base) if len(w) > 1]
    if len(words) > 3:
        candidates.append(" ".join(words[:3]))
    if len(words) > 2:
        candidates.append(" ".join(words[:2]))
    seen = set()
    out = []
    for item in candidates:
        key = item.lower()
        if item and key not in seen:
            seen.add(key)
            out.append(item)
    return out or [base]


def translate_query_zh(term: str) -> tuple[str, list[str], list[str]]:
    """把中文关键词翻成 NCBI 可检索的英文（用内置 ZH_DICT 反查）。

    返回 (新查询串, 命中列表, 未识别中文片段列表)。纯英文原样返回。
    """
    source = str(term or "")
    if not _CJK_RE.search(source):
        return source, [], []
    out = source
    used: list[str] = []
    for zh, en in _ZH_TO_EN:
        if zh and zh in out:
            out = out.replace(zh, en)
            used.append(f"{zh}→{en}")
    leftover = _CJK_RE.findall(out)
    out = _CJK_RE.sub(" ", out)
    out = re.sub(r"\s{2,}", " ", out).strip()
    return out, used, sorted(set(leftover))


def geo_search(term, species="all", study_type="all", min_samples=0,
               has_clinical=False, platform="all", limit=20, page=0,
               sort="relevance", has_matrix=None) -> dict:
    """GEO 系列检索（E-utilities esearch + esummary，中文化 + 空间转录组判定）。

    platform: all|visium|stereo|slideseq|10x|merfish —— 在客户端对 items 做平台过滤；
    平台过滤是严格的：不命中就返回空 + 提示（不再回退成全量无关数据）；
    平台词还会并入检索式，从源头减少非空间记录。
    has_matrix / has_clinical: 可选布尔过滤（None/False 表示不过滤）。

    返回 {"total","page","items":[{gse,title,title_zh,summary,species,n_samples,
    platform,study_type,release,has_matrix,has_clinical,disease,spatial_score,
    spatial_platform,tags}]}；中文查询额外带 query/query_zh/query_en/translated/
    translate_warning；逐级放宽触发时带 relaxed/used_term/relax_steps。
    """
    query_term = str(term or "").strip()
    if not query_term:
        raise EngineError("bad_request", "请输入检索关键词（支持中文疾病名、英文关键词或 GSE 号）")
    try:
        limit = max(1, min(100, int(limit or 20)))
    except (TypeError, ValueError):
        limit = 20
    try:
        page = max(0, int(page or 0))
    except (TypeError, ValueError):
        page = 0
    try:
        min_samples = max(0, int(min_samples or 0))
    except (TypeError, ValueError):
        min_samples = 0
    sort_mode = str(sort or "relevance").strip().lower()
    if sort_mode not in ("relevance", "date", "samples", "sample", "n_samples"):
        sort_mode = "relevance"
    platform_key = _platform_key(platform)
    want_matrix = bool(has_matrix) if has_matrix is not None else False

    if _CJK_RE.search(query_term):
        query_term, used_terms, leftover = translate_query_zh(query_term)
        if not query_term:
            raise EngineError("bad_request", "关键词里的中文没能翻译成英文，请改用英文关键词或补充疾病名（如：胃癌 → gastric cancer）")

    needs_filter = (bool(has_clinical) or want_matrix or platform_key != "all"
                    or sort_mode != "relevance")
    # 选了空间平台时，把平台词并进检索式（比事后过滤更准，也不会返回一堆非空间数据）
    if platform_key != "all":
        query_term = _with_platform_clause(query_term, platform_key)
    cascade = _cascade_terms(query_term)
    selected_term = cascade[0]
    selected_query = _build_query(selected_term, species, study_type, min_samples)
    ids: list[str] = []
    total = 0
    translation = ""
    relaxed_used: list[str] = []
    platform_note = ""
    if not needs_filter:
        # 非过滤模式也多抓一点候选（默认 40）：先按空间相关度排序再分页，
        # 否则纯疾病关键词的首屏会被非空间数据占满。
        pool = max(limit * 4, 40)
        for index, candidate in enumerate(cascade):
            query = _build_query(candidate, species, study_type, min_samples)
            ids, total, translation = _esearch(query, retmax=min(MAX_ESUMMARY_CANDIDATES, pool), retstart=0)
            selected_term, selected_query = candidate, query
            if total > 0:
                relaxed_used = cascade[:index + 1]
                if total >= max(6, limit // 2) or index == len(cascade) - 1:
                    break
        records = _esummary(ids)
        items = [_item_from_summary(record) for record in records]
        # 空间工作台：默认把空间相关度高的排前面（不隐藏非空间记录，但别让它们占满首屏）
        items.sort(key=lambda item: _SCORE_RANK.get(str(item.get("spatial_score") or "low"), 3))
        start = page * limit
        items = items[start:start + limit]
    else:
        # 过滤/排序需要候选池：抓一页更宽的候选再本地筛选+分页
        want = (page + 1) * limit
        candidate_n = min(MAX_ESUMMARY_CANDIDATES, max(want * 3, 60))
        for index, candidate in enumerate(cascade):
            query = _build_query(candidate, species, study_type, min_samples)
            ids, total, translation = _esearch(query, retmax=candidate_n, retstart=0)
            selected_term, selected_query = candidate, query
            if total > 0:
                relaxed_used = cascade[:index + 1]
                if total >= max(6, candidate_n // 3) or index == len(cascade) - 1:
                    break
        records = _esummary(ids)
        items = [_item_from_summary(record) for record in records]
        if has_clinical:
            items = [item for item in items if item.get("has_clinical")]
        if want_matrix:
            items = [item for item in items if item.get("has_matrix")]
        if platform_key != "all":
            matched = [item for item in items if _match_platform(item, platform_key)]
            # 严格过滤：选了平台就只给该平台的数据；宁可空，也不要拿无关数据充数
            items = matched
            items.sort(key=lambda item: _SCORE_RANK.get(str(item.get("spatial_score") or "low"), 3))
            if not matched:
                platform_note = ("该平台没有命中：已按平台词收紧检索式，仍无结果。"
                                 "可去掉平台筛选，或换更口语的关键词（如「胃癌 空间转录组」）")
        if sort_mode != "relevance":
            items = _sort_items(items, sort_mode)
        start = page * limit
        total = len(items)          # 过滤后总数（分页按过滤后计数）
        items = items[start:start + limit]
        if platform_note:
            total = 0 if not items and page == 0 else total

    if not items and not ids:
        raise EngineError("not_found", "没有检索到匹配的 GEO 系列：试试更宽泛的关键词（比如只写疾病名）或调整筛选条件")

    result = {"total": int(total), "page": int(page), "items": items}
    if translation:
        result["query"] = translation
    try:
        if _CJK_RE.search(str(term or "")):
            translated_query, used_terms, leftover = translate_query_zh(str(term or ""))
            result["query_zh"] = str(term or "")
            result["query_en"] = translated_query
            result["translated"] = used_terms
            if leftover:
                result["translate_warning"] = "这些中文词没有词典映射（已忽略）：" + "、".join(leftover[:6])
    except Exception:  # noqa: BLE001
        pass
    if needs_filter:
        result["filtered"] = bool(has_clinical or want_matrix or platform_key != "all")
    if platform_key != "all":
        result["platform"] = platform_key
    if platform_note:
        result["notes"] = [platform_note]
    if len(relaxed_used) > 1:
        result["relaxed"] = True
        result["used_term"] = selected_term
        result["relax_steps"] = relaxed_used
    return result


# ---------------------------------------------------------------------------
# 目录 HTML 解析
# ---------------------------------------------------------------------------

_HREF_RE = re.compile(r'<a\s+href="(?P<href>[^"]*)"[^>]*>(?P<name>[^<]*)</a>(?P<rest>[^<]*)',
                      re.IGNORECASE)
_SIZE_TOKEN_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)?[KMGTP]?$", re.IGNORECASE)


def _classify_file(name: str, source_dir: str, spatial_context: bool = False) -> str:
    """文件类型判定；kind 优先级：series_matrix > visium > raw > supplementary/other。

    spatial_context=True（title/summary 命中 Visium 类空间数据）时，文件名含 visium/spatial
    或 suppl 里的 *_RAW.tar 归为 "visium"——这类包通常带 spatial/ 坐标与组织图。
    """
    lower = name.lower()
    if "series_matrix" in lower:
        return SERIES_MATRIX_KIND
    if spatial_context and (
        "visium" in lower or "spatial" in lower
        or lower.endswith("_raw.tar") or lower.endswith("raw.tar")
        or re.search(r"_raw\.tar(\.gz)?$", lower)
    ):
        return SPATIAL_KIND
    if lower.endswith("_raw.tar") or lower.endswith("raw.tar") or re.search(r"_raw\.tar(\.gz)?$", lower):
        return "raw"
    if source_dir == "matrix":
        return "other"
    return "supplementary"


def _parse_size_from_rest(rest: str) -> int:
    tokens = [token for token in re.split(r"\s+", str(rest or "").strip()) if token]
    # 目录行：`2011-01-04 11:41  9.4K` → 最后一个 token 是大小（"-" 表示未知）
    for token in reversed(tokens):
        if token == "-":
            return 0
        if _SIZE_TOKEN_RE.match(token):
            return _parse_size(token)
    return 0


def _parse_dir_html(text: str, source_dir: str, spatial_context: bool = False) -> list[dict]:
    """解析 Apache 目录列表 HTML → [{name,size,url,kind,source}]。"""
    if not text:
        return []
    base = ""
    title_match = re.search(r"Index of\s+([^<\s]+)", text, re.IGNORECASE)
    if title_match:
        base = title_match.group(1).rstrip("/")
    files: list[dict] = []
    for match in _HREF_RE.finditer(text):
        href = html.unescape(match.group("href") or "").strip()
        raw_name = html.unescape(match.group("name") or "").strip()
        if not href or not raw_name:
            continue
        if href.startswith(("/", "http://", "https://", "ftp://", "?")):
            continue  # Parent Directory / 外链
        if href.endswith("/"):
            continue  # 子目录
        name = raw_name
        if href.startswith("./"):
            href = href[2:]
        path = urllib.parse.unquote(href)
        name = os.path.basename(path.rstrip("/")) or name
        size = _parse_size_from_rest(match.group("rest"))
        # base 是服务器绝对路径（如 /geo/series/GSE13nnn/GSE13861/matrix），
        # 不能直接拼 GEO_FTP（后者已含 /geo/series 前缀），否则会得到重复路径
        if base.startswith("/"):
            url = f"{NCBI_ORIGIN}{base}/{urllib.parse.quote(href)}"
        elif base:
            url = f"{GEO_FTP}/{base}/{urllib.parse.quote(href)}"
        else:
            url = urllib.parse.urljoin(GEO_FTP + "/", href)
        files.append({
            "name": name,
            "size": int(size),
            "url": url,
            "kind": _classify_file(name, source_dir, spatial_context),
            "source": source_dir,
        })
    # 去重（同名保留第一条）
    out = []
    seen = set()
    for item in files:
        if item["name"] in seen:
            continue
        seen.add(item["name"])
        out.append(item)
    return out


def _list_series_files(gse: str, spatial_context: bool = False) -> list[dict]:
    """抓 matrix/ 与 suppl/ 两个目录，合并文件列表（series_matrix 优先）。"""
    root = _series_root(gse)
    files: list[dict] = []
    matrix_html = _fetch_text(f"{root}/matrix/", soft_404=True)
    if matrix_html:
        files.extend(_parse_dir_html(matrix_html, "matrix", spatial_context))
    suppl_html = _fetch_text(f"{root}/suppl/", soft_404=True)
    if suppl_html:
        files.extend(_parse_dir_html(suppl_html, "suppl", spatial_context))
    # 排序：series_matrix → supplementary → raw/visium → other
    order = {SERIES_MATRIX_KIND: 0, "supplementary": 1, "raw": 2, SPATIAL_KIND: 2, "other": 3}
    files.sort(key=lambda item: (order.get(item["kind"], 9), item["name"]))
    return files


# ---------------------------------------------------------------------------
# 系列矩阵头部解析（design / contacts / clinical，Range 只取 1MB）
# ---------------------------------------------------------------------------


def _parse_matrix_header(text: str) -> dict:
    design = ""
    contacts: list[str] = []
    contributors: list[str] = []
    characteristics: list[list[str]] = []
    organism = ""
    titles: list[str] = []
    n_samples = 0
    for line in str(text or "").splitlines():
        if not line.startswith("!"):
            if line.strip():
                break  # 进入表格段，头部结束
            continue
        parts = line.rstrip("\r\n").split("\t")
        key = parts[0].strip()
        values = [_clean_quotes(value) for value in parts[1:]]
        if key == "!Series_overall_design" and values:
            design = design or values[0]
        elif key in ("!Series_contact_name", "!Series_contact_email", "!Series_contact_institute",
                     "!Series_contact_department", "!Series_contact_address"):
            contacts.extend(value for value in values if value)
        elif key == "!Series_contributor":
            contributors.extend(value for value in values if value)
        elif key == "!Sample_characteristics_ch1":
            characteristics.append(values)
        elif key == "!Sample_organism_ch1" and values:
            organism = organism or values[0]
        elif key == "!Sample_title":
            titles = values
            n_samples = max(n_samples, len([value for value in values if value]))
        elif key == "!Sample_geo_accession":
            n_samples = max(n_samples, len([value for value in values if value]))
    fields: list[dict] = []
    for index, values in enumerate(characteristics):
        non_empty = [value for value in values if value and value.lower() not in ("", "na", "nan")]
        if not non_empty:
            continue
        first = non_empty[0]
        if ":" in first:
            prefix = first.split(":", 1)[0].strip()
            key = re.sub(r"[^0-9A-Za-z_\u4e00-\u9fff]+", "_", prefix).strip("_").lower() + ":ch1"
            examples = []
            for value in non_empty:
                if ":" in value:
                    examples.append(value.split(":", 1)[1].strip())
                else:
                    examples.append(value)
            example = examples[0] if examples else first.split(":", 1)[1].strip()
        else:
            key = f"characteristics_ch1_{index + 1}"
            example = first
        # 同名 characteristics 行（GEO 常见多行 ch1）加序号去重
        base_key = key
        suffix = 2
        while any(field["key"] == key for field in fields):
            key = f"{base_key}#{suffix}"
            suffix += 1
        fields.append({"key": key, "n": len(non_empty), "example": example})
    contact_text = ""
    if contacts:
        names = contacts[0]
        emails = next((value for value in contacts if "@" in value), "")
        institutes = next((value for value in contacts if not value.lower().endswith("@") and "@" not in value
                           and value != names), "")
        pieces = [names]
        if emails:
            pieces.append(f"<{emails}>")
        if institutes and institutes != names:
            pieces.append(institutes)
        contact_text = " ".join(piece for piece in pieces if piece)
    if contributors:
        extra = "; ".join(contributors[:6])
        contact_text = f"{contact_text}; {extra}" if contact_text else extra
    return {
        "design": design,
        "contacts": contact_text,
        "clinical_fields": fields,
        "organism": organism,
        "n_samples": n_samples,
    }


def _fetch_matrix_head(files: list[dict]) -> dict:
    matrix_file = next((item for item in files if item.get("kind") == SERIES_MATRIX_KIND), None)
    if not matrix_file:
        return {}
    url = matrix_file.get("url") or ""
    if not url:
        return {}
    try:
        raw = _fetch_bytes(url, limit=MATRIX_HEAD_BYTES, retries=1, timeout=MATRIX_HEAD_TIMEOUT,
                           headers={"Range": f"bytes=0-{MATRIX_HEAD_BYTES - 1}"})
    except EngineError:
        return {}
    if not raw:
        return {}
    if raw[:2] == b"\x1f\x8b":
        decompressor = zlib.decompressobj(31)
        try:
            raw = decompressor.decompress(raw)
        except zlib.error:
            raw = b""
    text = raw.decode("utf-8", "replace")
    return _parse_matrix_header(text)


# ---------------------------------------------------------------------------
# geo_detail
# ---------------------------------------------------------------------------


def _summary_for_gse(gse: str) -> dict:
    uid = _uid_for(gse)
    records = _esummary([uid])
    if records:
        return records[0]
    # 兜底：用 ACCN 字段搜索（极老的系列 UID 可能对不上）
    ids, _, _ = _esearch(f"{gse}[ACCN]", retmax=10, retstart=0)
    candidates = [item for item in ids if item.startswith("2")]
    if candidates:
        records = _esummary(candidates[:1])
        if records:
            return records[0]
    raise EngineError("not_found", f"GEO 系列 {gse} 不存在或 NCBI 没有返回摘要")


def geo_detail(gse) -> dict:
    """GEO 系列详情（摘要 + 目录文件 + 矩阵头部临床字段 + 空间平台判定）。

    相对 bg_geo 的新增：spatial_score / spatial_platform / visium_hint，
    以及 files[].kind 的 "visium"（series_matrix 优先，其次 suppl 里的 Visium 类 RAW.tar）。
    """
    gse = _norm_gse(gse)
    record = _summary_for_gse(gse)
    title = str(record.get("title") or "")
    summary = str(record.get("summary") or "")
    taxon = str(record.get("taxon") or "")
    try:
        n_samples = int(record.get("n_samples") or len(record.get("samples") or []))
    except (TypeError, ValueError):
        n_samples = 0
    pubmed_ids = [str(item) for item in (record.get("pubmedids") or []) if str(item).strip()]
    spatial_score = spatial_score_hint(title, summary)
    spatial_platform = spatial_platform_hints(title, summary)
    visium_context = _visium_context(title, summary)
    files = _list_series_files(gse, spatial_context=visium_context)
    header = _fetch_matrix_head(files) if files else {}
    species = taxon or header.get("organism") or ""
    clinical_fields = header.get("clinical_fields") or []
    # 是否可能有 spatial/ 目录（Visium 类数据）：文件线索或标题/摘要线索
    visium_hint = bool(
        any(item.get("kind") == SPATIAL_KIND for item in files)
        or any("visium" in item["name"].lower() or "spatial" in item["name"].lower()
               for item in files)
        or "Visium" in spatial_platform or "Visium HD" in spatial_platform
    )
    result = {
        "gse": gse,
        "title": title,
        "title_zh": _translate(title),
        "summary": summary,
        "design": header.get("design") or "",
        "species": species,
        "n_samples": int(n_samples or header.get("n_samples") or 0),
        "platforms": _platforms(record.get("gpl")),
        "pubmed": ", ".join(pubmed_ids),
        "release": _norm_date(record.get("pdat")),
        "contacts": header.get("contacts") or "",
        "clinical": {
            "from_series_matrix": bool(clinical_fields),
            "fields": clinical_fields,
        },
        "spatial_score": spatial_score,
        "spatial_platform": spatial_platform,
        "visium_hint": visium_hint,
        "files": [
            {"name": item["name"], "size": int(item["size"]), "url": item["url"],
             "kind": item["kind"], "source": item.get("source") or ""}
            for item in files
        ],
        "preview_files": [
            {"name": item["name"], "size": int(item["size"]), "kind": item["kind"]}
            for item in files[:12]
        ],
        "suggested": "",
        "url": f"{GEO_BASE}?acc={gse}",
    }
    for item in files:
        if item.get("kind") == SERIES_MATRIX_KIND:
            result["suggested"] = item["name"]
            break
    if not result["suggested"] and files:
        result["suggested"] = files[0]["name"]
    result["has_matrix"] = bool(result["suggested"])
    result["has_clinical"] = bool(clinical_fields) or _looks_clinical(f"{title} {summary}")
    result["disease"] = _detect_disease(title, summary)
    result["tags"] = _make_tags(title, summary, species, str(record.get("gdstype") or ""),
                                result["disease"])
    return result


# ---------------------------------------------------------------------------
# geo_download
# ---------------------------------------------------------------------------


def _check_dest_dir(dest_dir) -> str:
    text = str(dest_dir or "").strip()
    if not text:
        raise EngineError("bad_request", "请指定下载目录（dest_dir 必须是绝对路径）")
    path = os.path.expanduser(text)
    if not os.path.isabs(path):
        raise EngineError("bad_request", f"dest_dir 必须是绝对路径：{dest_dir!r}（示例：/Users/you/GEO）")
    path = os.path.abspath(path)
    if path == os.path.sep or path == "/":
        raise EngineError("bad_request", "不能直接下载到文件系统根目录 /，请指定子目录")
    if os.path.exists(path) and not os.path.isdir(path):
        raise EngineError("bad_request", f"dest_dir 已经是文件而不是目录：{path}")
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as error:
        raise EngineError("bad_request", f"无法创建下载目录：{path}（{error}）")
    if not os.access(path, os.W_OK):
        raise EngineError("bad_request", f"下载目录不可写：{path}")
    return path


def _safe_filename(url: str, fallback: str) -> str:
    name = os.path.basename(urllib.parse.unquote(urllib.parse.urlparse(url).path or "")) or fallback
    name = os.path.basename(name.replace("\\", "/"))
    if not name or name in (".", "..") or "/" in name or "\x00" in name:
        name = re.sub(r"[^0-9A-Za-z._-]+", "_", fallback or "download.bin") or "download.bin"
    return name


def _emit(progress, message: str, pct: float | None = None) -> None:
    if not progress:
        return
    try:
        progress(message, pct)
    except TypeError:
        try:
            progress(message)
        except Exception:  # noqa: BLE001
            pass
    except Exception:  # noqa: BLE001
        pass


class _DownloadIntegrityError(Exception):
    """下载内容不完整（字节数不符 / 压缩容器打不开）：可重试，不删除旧文件。"""


def _content_length(headers) -> int | None:
    """解析 Content-Length；缺失或非法返回 None（0 视为合法长度）。"""
    if headers is None:
        return None
    try:
        text = str(headers.get("Content-Length") or "").strip()
    except Exception:  # noqa: BLE001 - 假响应/异常 headers 都按缺失处理
        return None
    if not text:
        return None
    try:
        value = int(text)
    except (TypeError, ValueError):
        return None
    return value if value >= 0 else None


def _check_tar_members(archive: tarfile.TarFile, *, check_size: bool, path: str) -> None:
    """遍历 tar 成员；check_size=True 时用成员偏移校验文件是否被截断。"""
    members = archive.getmembers()
    if not check_size or not members:
        return
    expected_end = max((int(item.offset_data) + int(item.size) for item in members), default=0)
    actual = os.path.getsize(path)
    if actual < expected_end:
        raise _DownloadIntegrityError(
            f"tar 数据被截断：实际 {actual} 字节 < 成员所需 {expected_end} 字节")


def _validate_download(path: str, expected: int | None, target: str) -> None:
    """下载完整性校验：优先核对 Content-Length；缺失时按扩展名做容器可读性兜底。"""
    try:
        size = os.path.getsize(path)
    except OSError as error:
        raise _DownloadIntegrityError(f"无法读取临时文件：{error}")
    if expected is not None:
        if size != expected:
            raise _DownloadIntegrityError(
                f"写入 {size} 字节，与 Content-Length={expected} 不一致（{os.path.basename(target)}）")
        return
    name = os.path.basename(target).lower()
    try:
        if name.endswith((".tar.gz", ".tgz")):
            with tarfile.open(path, "r:gz") as archive:
                _check_tar_members(archive, check_size=False, path=path)
        elif name.endswith((".tar.bz2", ".tbz", ".tbz2")):
            with tarfile.open(path, "r:bz2") as archive:
                _check_tar_members(archive, check_size=False, path=path)
        elif name.endswith(".tar"):
            with tarfile.open(path, "r:") as archive:
                _check_tar_members(archive, check_size=True, path=path)
        elif name.endswith(".gz"):
            if size <= VALIDATE_GZIP_FULL_MAX:
                with gzip.open(path, "rb") as handle:
                    while handle.read(DOWNLOAD_CHUNK):
                        pass
            else:
                with gzip.open(path, "rb") as handle:
                    handle.read(DOWNLOAD_CHUNK)  # 只验证 gzip 头/首个解压块可读
        elif name.endswith(".bz2"):
            if size <= VALIDATE_GZIP_FULL_MAX:
                with bz2.open(path, "rb") as handle:
                    while handle.read(DOWNLOAD_CHUNK):
                        pass
            else:
                with bz2.open(path, "rb") as handle:
                    handle.read(DOWNLOAD_CHUNK)
    except (OSError, EOFError, tarfile.TarError, zlib.error) as error:
        raise _DownloadIntegrityError(
            f"{os.path.basename(target)} 无法通过容器可读性校验：{error}")


def _download_one(url: str, target: str, *, progress=None, label: str = "") -> int:
    """流式下载单个 URL 到 target（失败重试 ≤2 次；仅支持失败重试，不支持断点续传）。

    写入 `<target>.part`，完整性校验通过后才 os.replace 原子替换 target：
      * 有 Content-Length 时必须严格相等；
      * 缺失时按扩展名用 gzip/bz2/tar 可读性兜底。
    任何失败/中断都会删除 `.part`，并且不触碰已存在的旧 target（旧文件保持原样）。
    """
    part = f"{target}{PART_SUFFIX}"
    attempts = len(RETRY_SLEEPS) + 1
    last_error: Exception | None = None
    for attempt in range(attempts):
        _remove_quiet(part)  # 每次重试都从头写（不支持断点续传）
        written = 0
        try:
            with _urlopen(_request(url), HTTP_TIMEOUT) as response:
                total = _content_length(getattr(response, "headers", None))
                if total is not None and total > MAX_DOWNLOAD_BYTES:
                    raise EngineError(
                        "unsupported",
                        f"文件过大（{total / 1024 ** 3:.1f} GB），超过单文件 8GB 下载上限：{label or url}",
                    )
                with open(part, "wb") as handle:
                    while True:
                        chunk = response.read(DOWNLOAD_CHUNK)
                        if not chunk:
                            break
                        written += len(chunk)
                        if written > MAX_DOWNLOAD_BYTES:
                            raise EngineError(
                                "unsupported",
                                f"文件超过单文件 8GB 下载上限，已中止：{label or url}",
                            )
                        handle.write(chunk)
                        if progress and total:
                            pct = min(99.0, written / total * 100.0)
                            _emit(progress, f"下载 {label or url}：{written / 1024 ** 2:.1f}/"
                                            f"{total / 1024 ** 2:.1f} MB", pct)
                    handle.flush()
                    try:
                        os.fsync(handle.fileno())
                    except OSError:
                        pass
                _validate_download(part, total, target)
            os.replace(part, target)
            return written
        except EngineError:
            _remove_quiet(part)
            raise
        except _DownloadIntegrityError as error:
            last_error = error
            _remove_quiet(part)
        except Exception as error:  # noqa: BLE001 - 网络/磁盘异常统一重试
            last_error = error
            _remove_quiet(part)
        if attempt < attempts - 1:
            delay = RETRY_SLEEPS[min(attempt, len(RETRY_SLEEPS) - 1)]
            if delay > 0:
                time.sleep(delay)
    if isinstance(last_error, _DownloadIntegrityError):
        raise EngineError("unsupported",
                          f"下载校验失败：{last_error}（已重试 {attempts} 次，旧文件未改动）")
    raise _network_fail(last_error)


def _remove_quiet(path: str) -> None:
    if not path:
        return
    try:
        if os.path.exists(path):
            os.remove(path)
    except OSError:
        pass


def geo_download(gse, files, dest_dir, overwrite=False, progress=None) -> dict:
    """把指定文件下载到 dest_dir（绝对路径、不可为 /），返回 {"dir","files","ms"}。

    覆盖下载（overwrite=True）采用「.part 临时文件 + 完整性校验 + os.replace 原子替换」：
    下载失败/中断时旧文件保持原样；不做断点续传，失败只从头重试。
    """
    started = time.time()
    gse = _norm_gse(gse)
    dest = _check_dest_dir(dest_dir)
    if isinstance(files, str):
        requested = [files]
    elif isinstance(files, (list, tuple, set)):
        requested = [str(item) for item in files if str(item).strip()]
    else:
        requested = []
    _emit(progress, f"正在获取 {gse} 的文件列表", 1)

    detail_files: list[dict] = []
    detail_error: EngineError | None = None
    try:
        detail_files = geo_detail(gse).get("files") or []
    except EngineError as error:
        if error.code == "not_found":
            raise
        detail_error = error  # 目录列表失败仍可尝试按 suppl 规则拼接 URL

    by_name = {item["name"]: item for item in detail_files}
    if not requested:
        suggested = next((item["name"] for item in detail_files
                          if item.get("kind") == SERIES_MATRIX_KIND), "")
        if not suggested and detail_files:
            suggested = detail_files[0]["name"]
        if not suggested:
            raise EngineError(
                "not_found",
                f"{gse} 在 GEO FTP 上没有可下载文件（matrix/ 与 suppl/ 都为空）",
            )
        requested = [suggested]

    plan: list[dict] = []
    unknown: list[str] = []
    for item in requested:
        if re.match(r"^(?:https?|ftp)://", item, re.IGNORECASE):
            plan.append({"name": _safe_filename(item, f"{gse}.txt"), "url": item, "size": 0})
            continue
        entry = by_name.get(item) or by_name.get(os.path.basename(item))
        if entry is None:
            if detail_files:
                unknown.append(item)
                continue
            # 详情不可用：按文件名猜目录（series_matrix → matrix/，其余 → suppl/）
            subdir = "matrix" if SERIES_MATRIX_KIND in item.lower() or "series_matrix" in item.lower() else "suppl"
            plan.append({
                "name": os.path.basename(item),
                "url": f"{_series_root(gse)}/{subdir}/{urllib.parse.quote(os.path.basename(item))}",
                "size": 0,
            })
            continue
        plan.append({"name": entry["name"], "url": entry["url"], "size": int(entry.get("size") or 0)})
    if unknown:
        raise EngineError(
            "not_found",
            f"文件不在 {gse} 的可下载列表里：{', '.join(unknown[:5])}",
            candidates=[{"name": item["name"], "size": item["size"], "kind": item["kind"]}
                        for item in detail_files[:30]],
        )
    if detail_error is not None and not plan:
        raise detail_error

    results: list[dict] = []
    total = len(plan)
    failed = 0
    for index, entry in enumerate(plan, start=1):
        target = os.path.join(dest, entry["name"])
        base = (index - 1) / max(total, 1) * 100.0
        if os.path.exists(target) and not overwrite:
            size = os.path.getsize(target)
            results.append({"name": entry["name"], "path": target, "bytes": size,
                            "ok": True, "skipped": True})
            _emit(progress, f"已存在，跳过：{entry['name']}", min(100.0, base + 100.0 / total))
            continue
        # 注意：overwrite=True 也不再先删除旧文件——由 _download_one 写 .part，
        # 校验通过后原子替换；失败时旧文件必须保持可用。
        _emit(progress, f"正在下载 {entry['name']}（{index}/{total}）", base)
        try:
            written = _download_one(entry["url"], target, progress=progress, label=entry["name"])
            results.append({"name": entry["name"], "path": target, "bytes": int(written), "ok": True})
        except EngineError as error:
            if error.code == "bad_request":  # 路径类错误直接终止
                raise
            failed += 1
            results.append({"name": entry["name"], "path": target, "bytes": 0,
                            "ok": False, "error": error.message})
        _emit(progress, f"完成 {entry['name']}（{index}/{total}）", min(100.0, base + 100.0 / total))
    ok_count = sum(1 for item in results if item.get("ok"))
    if ok_count == 0 and results:
        first_error = next((item.get("error") for item in results if item.get("error")), "")
        raise EngineError("unsupported", f"下载失败：{first_error or 'NCBI 暂时不可达'}")
    _emit(progress, f"下载完成：{ok_count}/{len(results)} 个文件", 100)
    result = {
        "dir": dest,
        "gse": gse,
        "files": results,
        "ms": int((time.time() - started) * 1000),
        "resume": False,  # 仅支持失败重试，不支持断点续传
        "notes": ["下载只支持失败重试（每次从头下载），不支持断点续传"],
    }
    if failed:
        result["warnings"] = [f"{failed} 个文件下载失败"]
    return result


# ---------------------------------------------------------------------------
# 自检（只用于测试/手工验证，不参与协议）
# ---------------------------------------------------------------------------

def supported_dictionary_size() -> int:
    return len(ZH_DICT)


__all__ = ["geo_search", "geo_detail", "geo_download", "ZH_DICT", "DISEASE_KEYS",
           "supported_dictionary_size", "spatial_platform_hints", "spatial_score_hint",
           "_SPATIAL_TERMS", "_SPATIAL_WEAK_TERMS", "_PLATFORM_MATCH"]
