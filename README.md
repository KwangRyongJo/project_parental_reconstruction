###########parental haplotype reconstruction########

###script_Phasing_PCA_HC_Fst

############phasing##########

# 마스터 파이프라인 실행
bash 00-run_all.sh

#output folders/files
haploblock_out/
haplotypes_out/
lists/
logs/
out_parents/
reports/결과파일 있음
input_data.plink.bed
input_data.plink.bim
input_data.plink.fam
input_data.plink.log
input_data.plink.map
input_data.plink.nosex
input_data.plink.ped
qc_check.frq
qc_check.imiss
qc_check.lmiss
qc_check.log
qc_check.nosex

###########PCA_HC_Fst_LD#################

setwd("D:/2026/scripts/parent_infer/PCA_HC")

# ---- 패키지 없이 base R로도 가능하지만, 깔끔한 CSV 입출력을 위해 readr 사용(있으면) ----
use_readr <- requireNamespace("readr", quietly = TRUE)

# 입력/출력 파일 이름
infile  <- "input.csv"
outfile <- "output.txt"

# ---------- 데이터 읽기 ----------
if (use_readr) {
  df <- readr::read_csv(infile, show_col_types = FALSE)
} else {
  df <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
}

# 필수 열 확인
req_cols <- c("ID", "chrom", "position", "ref")
missing <- setdiff(req_cols, names(df))
if (length(missing) > 0) {
  stop(paste("다음 필수 열이 없습니다:", paste(missing, collapse = ", ")))
}

# 샘플(유전자형) 열 자동 식별
sample_cols <- setdiff(names(df), req_cols)
if (length(sample_cols) == 0) {
  stop("샘플 유전자형 열(P01, P02, ...)이 없습니다.")
}

# ---------- chr 번호 추출 (Ca_59Chr01 ~ Ca_59Chr12 → 1 ~ 12) ----------
# chrom 끝 두 자리 숫자를 정수로 변환
chr_str <- gsub(".*([0-9]{2})$", "\\1", df$chrom)
chr_num <- suppressWarnings(as.integer(chr_str))
if (any(is.na(chr_num))) {
  stop("chrom에서 말미 숫자(예: Ca_59Chr01) 추출 실패한 행이 있습니다.")
}

# ---------- ALT(대립유전자) 계산: 단일 alt 가정 ----------
# 각 행의 샘플 유전자형들(P01~)에서 ref가 아닌 염기를 모아 사전순으로 하나 선택
find_alt_single <- function(ref_base, geno_row_vec) {
  # 유전자형을 글자 단위로 분해: 예 "CA" → c("C","A")
  chars <- unlist(strsplit(paste0(geno_row_vec, collapse = ""), split = ""))
  chars <- toupper(chars)
  # 염기만 필터
  bases <- chars[chars %in% c("A","C","G","T")]
  # ref 제외
  alts <- setdiff(unique(bases), toupper(ref_base))
  if (length(alts) == 0) NA_character_ else sort(alts)[1]
}

alt_vec <- vapply(
  seq_len(nrow(df)),
  function(i) {
    find_alt_single(df$ref[i], as.character(unlist(df[i, sample_cols])))
  },
  FUN.VALUE = character(1)
)

# ---------- 출력 데이터 구성 ----------
# position을 가능한 정수로 변환(문자 섞이면 NA 가능)
pos_int <- suppressWarnings(as.integer(df$position))

marker  <- paste0("Chr", chr_num, "_", pos_int)
alleles <- ifelse(is.na(alt_vec),
                  paste0(toupper(df$ref), "/", toupper(df$ref)),  # alt 없음 → ref/ref로 표시
                  paste0(toupper(df$ref), "/", alt_vec))

out_df <- data.frame(
  marker  = marker,
  alleles = alleles,
  chr     = chr_num,
  pos     = pos_int,
  df[, sample_cols, drop = FALSE],
  check.names = FALSE
)

# ---------- 저장 ----------
if (use_readr) {
  # readr 사용: 탭으로 구분된 txt 저장 (UTF-8)
  readr::write_tsv(out_df, outfile, na = "")
} else {
  # base R 사용: 탭으로 구분된 txt 저장 (엑셀 호환, UTF-8 BOM)
  write.table(
    out_df,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    fileEncoding = "UTF-8-BOM"
  )
}

message("저장 완료: ", normalizePath(outfile))

########################## 
setwd("D:/2026/scripts/parent_infer/PCA_HC")

# 데이터 불러오기
data <- read.table("output.txt", header = TRUE, sep = "\t")
#long time!!!
# 유전자형 매핑 함수 정의
map_genotype <- function(genotype, alleles) {
  if (is.na(alleles) || alleles == "") {
    return(NA) # alleles가 NA이거나 빈 문자열일 경우 NA 반환
  }
  
  # alleles가 팩터일 경우 문자열로 변환
  if (is.factor(alleles)) {
    alleles <- as.character(alleles)
  }
  
  ref <- unlist(strsplit(alleles, "/"))[1]
  alt <- unlist(strsplit(alleles, "/"))[2]
  
  if (is.na(genotype)) {
    return(NA)
  } else if (genotype == paste0(ref, ref)) {
    return(-1)
  } else if (genotype == paste0(alt, alt)) {
    return(1)
  } else if (genotype == paste0(ref, alt) || genotype == paste0(alt, ref)) {
    return(0)
  } else {
    return(NA) # 결측값 처리
  }
}

# 각 유전자형 데이터를 매핑
for (i in 1:nrow(data)) {
  alleles <- data$alleles[i]
  if (is.null(alleles) || alleles == "") {
    next # alleles가 NULL이거나 빈 문자열이면 다음 반복으로 넘어감
  }
  for (j in 5:ncol(data)) { # 첫 4개 컬럼은 메타데이터
    data[i, j] <- map_genotype(data[i, j], alleles)
  }
}

# 불필요한 열 제거
#data <- data[, !names(data) %in% c("alleles", "chr", "pos")]

# 결과 저장
write.table(data, "output_numerical.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################### 


#############################################
## Ca 유전형(.txt) + 지도(.csv) 전체 파이프라인
## - 결측률 필터링
## - 지도 기준 정렬(ChrX + 위치)
## - rrBLUP 보간(EM) & Kinship
## - PCA / Screeplot / PC1-2 산점도
## - 계층적 클러스터링 (Ward.D2)
## - 결과 저장(PDF/CSV)
#############################################

## 0) 환경 설정 ------------------------------------------------------
setwd("D:/2026/scripts/parent_infer/PCA_HC")
options(stringsAsFactors = FALSE)
set.seed(123)

## 1) 패키지 로드 ----------------------------------------------------
library(rrBLUP)
library(corrplot)   # 경고(작성 R 버전)는 무시 가능
# library(ggplot2)  # 필요시 사용

## 2) 데이터 로드 ----------------------------------------------------
# .txt 유전형: 탭 구분 + 첫 열(마커ID)을 행이름으로 사용
geno <- read.table("output_numerical.txt",
                   header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE)

# .csv 지도: 인덱스 숫자 컬럼이 있으므로 row.names=1을 쓰지 말고,
# 실제 마커ID가 들어있는 'Markers' 컬럼로 행이름을 지정
map  <- read.csv("Cs_map.csv", header = TRUE, check.names = FALSE)
stopifnot("Markers" %in% colnames(map))
map <- map[!is.na(map$Markers) & !duplicated(map$Markers), ]
rownames(map) <- as.character(map$Markers)

## 3) 메타 열 제거(alleles/chr/pos) 후 유전형 매트릭스 확보 ---------
meta_cols <- c("alleles", "chr", "pos")
geno_mat  <- geno[, setdiff(colnames(geno), meta_cols)]
cat("geno_mat dim =", dim(geno_mat), "\n")  # (마커수 x 샘플수) = (44325 x 205) 예상

## 4) 결측률 필터링 ---------------------------------------------------
# - 행=마커, 열=샘플 구조에서 결측률 계산
filter.fun <- function(geno, IM = 0.40, MM = 0.60) {
  inds_miss  <- colMeans(is.na(geno))                # 샘플 결측률(열 기준)
  marks_miss <- rowMeans(is.na(geno))                # 마커 결측률(행 기준)

  keep_inds  <- which(inds_miss  < IM)
  keep_marks <- which(marks_miss < MM)

  geno_f <- geno[keep_marks, keep_inds, drop = FALSE]
  return(geno_f)
}

geno.filtered <- filter.fun(geno_mat, IM = 0.40, MM = 0.60)
cat("geno.filtered dim =", dim(geno.filtered), "\n")

## (선택) 필터 후 결측률 요약 저장
ind_miss  <- colMeans(is.na(geno.filtered))
mark_miss <- rowMeans(is.na(geno.filtered))
write.csv(data.frame(sample = colnames(geno.filtered), missing_rate = ind_miss),
          "samples_missing_rate.csv", row.names = FALSE)
write.csv(data.frame(marker = rownames(geno.filtered), missing_rate = mark_miss),
          "markers_missing_rate.csv", row.names = FALSE)

## 5) 지도 기준 정렬(Chr/좌표) ---------------------------------------
# 지도 컬럼 확인: 'chrom'(예: 'Chr1'), 'loc'(좌표)
stopifnot(all(c("chrom", "loc") %in% colnames(map)))

# 공통 마커만 사용
common_marks <- intersect(rownames(geno.filtered), rownames(map))
cat("공통 마커 수 =", length(common_marks), "\n")

map_sub <- map[common_marks, , drop = FALSE]

# 정렬용 숫자 염색체 열 생성: 'Chr1' -> 1, 'Chr12' -> 12
map_sub$chrom_num <- suppressWarnings(as.integer(sub("^Chr", "", map_sub$chrom)))
map_sub$loc       <- suppressWarnings(as.numeric(map_sub$loc))

# NA 행 제거 후 (chrom_num, loc) 기준 정렬
map_sub <- map_sub[complete.cases(map_sub[, c("chrom_num", "loc")]), , drop = FALSE]
ord     <- order(map_sub$chrom_num, map_sub$loc)
map_sub <- map_sub[ord, , drop = FALSE]

# 유전형 행(=마커)을 지도 순서로 재배치
geno.filtered <- geno.filtered[rownames(map_sub), , drop = FALSE]
cat("정렬 후 geno.filtered dim =", dim(geno.filtered), "\n")

## 6) rrBLUP 보간 & Kinship ------------------------------------------
# rrBLUP::A.mat은 행=샘플, 열=마커가 일반적 → 전치
geno_rr <- t(geno.filtered)  # (샘플 x 마커)

# A.mat 시도: EM 대치 → 실패 시 평균대치로 폴백
imp_try <- try(
  A.mat(geno_rr, impute.method = "EM",
        return.imputed = TRUE, min.MAF = 0.05),
  silent = TRUE
)

if (inherits(imp_try, "try-error")) {
  message("EM 대치 실패 → 평균대치로 재시도합니다.")
  Imputation <- A.mat(geno_rr, impute.method = "mean",
                      return.imputed = TRUE, min.MAF = 0.05)
} else {
  Imputation <- imp_try
}

K.mat     <- Imputation$A          # 근친행렬 (샘플 x 샘플)
geno.gwas <- Imputation$imputed    # 보간 유전형 (샘플 x 마커)

cat("geno.gwas dim =", dim(geno.gwas), "\n")
cat("K.mat dim =", dim(K.mat), "\n")

# 일부 확인
if (nrow(geno.gwas) >= 5 && ncol(geno.gwas) >= 5) {
  print(geno.gwas[1:5, 1:5])
}
if (nrow(K.mat) >= 5 && ncol(K.mat) >= 5) {
  print(K.mat[1:5, 1:5])
}

# 결과 저장
#write.csv(geno.gwas, "geno_imputed.csv")
write.table(geno.gwas, "geno_imputed.tsv")

write.csv(K.mat,     "kinship_matrix.csv")

# Kinship 히트맵 (corrplot)
pdf("Figure1_kinship.pdf")
corrplot(K.mat, tl.cex = 0.5, method = "color", tl.pos = "n",
         main = "Kinship heatmap (rrBLUP A.mat)")
dev.off()

## 7) PCA / Screeplot / PC1-2 산점도 -------------------------------
# 데이터 중심화(center=TRUE). scale=FALSE(표준화 없음)
geno.scale <- scale(geno.gwas, center = TRUE, scale = FALSE)

# 차원 체크
if (any(dim(geno.scale) == 0)) {
  stop("PCA를 위한 데이터 차원이 0입니다. 이전 단계(지도 정렬/MAF/필터링)를 확인하세요.")
}

svdgeno <- svd(geno.scale)
PCA     <- geno.scale %*% svdgeno$v  # (샘플 x PC)

# Screeplot (상위 10개 PC)
pdf("Figure2_merged.pdf")
plot(round((svdgeno$d)^2 / sum((svdgeno$d)^2), d = 7)[1:10],
     type = "o", main = "Screeplot",
     xlab = "PCAs", ylab = "% variance")
dev.off()

# PC1/PC2 기여율
PCA1 <- 100 * round((svdgeno$d[1])^2 / sum((svdgeno$d)^2), d = 3)
PCA2 <- 100 * round((svdgeno$d[2])^2 / sum((svdgeno$d)^2), d = 3)
cat(sprintf("PCA1 = %.3f %% | PCA2 = %.3f %%\n", PCA1, PCA2))

# PC1 vs PC2 산점도
pdf("Figure3_merged.pdf")
plot(PCA[, 1], PCA[, 2],
     xlab = paste("Pcomp:", PCA1, "%", sep = ""),
     ylab = paste("Pcomp:", PCA2, "%", sep = ""),
     pch = 20, cex = 0.7, main = "PC1 vs PC2")
dev.off()



## PCA 결과(PCA 객체)가 준비된 상태에서 색상 지정
# 샘플명 추출
samples <- rownames(geno.gwas)

# 색상 벡터 생성: P로 시작하면 빨간색, 아니면 회색
colors <- ifelse(grepl("^P", samples), "red", "gray")

# PC1/PC2 기여율 계산
PCA1 <- 100 * round((svdgeno$d[1])^2 / sum((svdgeno$d)^2), d = 3)
PCA2 <- 100 * round((svdgeno$d[2])^2 / sum((svdgeno$d)^2), d = 3)

# PDF로 저장
pdf("Figure3_merged_colored.pdf")
plot(PCA[, 1], PCA[, 2],
     col = colors, pch = 20, cex = 0.8,
     xlab = paste("PC1:", PCA1, "%", sep = ""),
     ylab = paste("PC2:", PCA2, "%", sep = ""),
     main = "PCA Plot (P 시작 샘플 = 빨간색)")
legend("topright", legend = c("P 시작 샘플", "기타 샘플"),
       col = c("red", "gray"), pch = 20)
dev.off()

========================= 
## --- 전제: geno.gwas (행=샘플, 열=마커)가 존재 ---
stopifnot(is.matrix(geno.gwas) || is.data.frame(geno.gwas))

## 패키지
# install.packages(c("ggplot2","ggrepel","showtext","sysfonts","Cairo")) # 필요시
library(ggplot2)
library(ggrepel)
library(showtext)
library(sysfonts)
library(Cairo)

## 0) 한글 폰트(네모 방지)
font_path_candidates <- c(
  "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
  "/usr/share/fonts/opentype/noto/NotoSansCJKkr-Regular.otf",
  "/usr/share/fonts/truetype/noto/NotoSansCJKkr-Regular.otf",
  "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
  "C:/Windows/Fonts/malgun.ttf"
)
font_path <- font_path_candidates[file.exists(font_path_candidates)][1]
if (!is.na(font_path)) {
  font_add("korean", regular = font_path)
  showtext_auto(TRUE)
}

## 1) P코드 -> 표시이름 매핑
map_P_to_name <- c(
  "P01"="굿모닝",
  "P01_ParentA"="굿모닝_Parent1",
  "P01_ParentB"="굿모닝_Parent2",
  "P02"="단청",
  "P02_ParentA"="단청_Parent1",
  "P02_ParentB"="단청_Parent2",
  "P03"="설향",
  "P03_ParentA"="설향_Parent1",
  "P03_ParentB"="설향_Parent2",
  "P04"="봄봄",
  "P04_ParentA"="봄봄_Parent1",
  "P04_ParentB"="봄봄_Parent2",  
  "P05"="청춘",
  "P05_ParentA"="청춘_Parent1",
  "P05_ParentB"="청춘_Parent2",
  "P06"="청어람",
  "P06_ParentA"="청어람_Parent1",
  "P06_ParentB"="청어람_Parent2",
  "P07"="백리향",
  "P07_ParentA"="백리향_Parent1",
  "P07_ParentB"="백리향_Parent2",
  "P08"="천리향",
  "P08_ParentA"="천리향_Parent1",
  "P08_ParentB"="천리향_Parent2",
  "P09"="기라성",
  "P09_ParentA"="기라성_Parent1",
  "P09_ParentB"="기라성_Parent2",
  "P10"="스트롱",
  "P10_ParentA"="스트롱_Parent1",
  "P10_ParentB"="스트롱_Parent2",
  "P11"="베테랑",
  "P11_ParentA"="베테랑_Parent1",
  "P11_ParentB"="베테랑_Parent2",
  "P12"="하이탑",
  "P12_ParentA"="하이탑_Parent1",
  "P12_ParentB"="하이탑_Parent2",
  "P13"="썸머킹",
  "P13_ParentA"="썸머킹_Parent1",
  "P13_ParentB"="썸머킹_Parent2",
  "P14"="진한백다다기",
  "P14_ParentA"="진한백다다기_Parent1",
  "P14_ParentB"="진한백다다기_Parent2",
  "P15"="낙동청장",
  "P15_ParentA"="낙동청장_Parent1",
  "P15_ParentB"="낙동청장_Parent2",
  "P16"="신동청장",
  "P16_ParentA"="신동청장_Parent1",
  "P16_ParentB"="신동청장_Parent2",
  "P17"="강풍삼척",
  "P17_ParentA"="강풍삼척_Parent1",
  "P17_ParentB"="강풍삼척_Parent2",
  "P18"="장형낙합",
  "P18_ParentA"="장형낙합_Parent1",
  "P18_ParentB"="장형낙합_Parent2",
  "P19"="신세대",
  "P19_ParentA"="신세대_Parent1",
  "P19_ParentB"="신세대_Parent2",
  "P20"="새로미",
  "P20_ParentA"="새로미_Parent1",
  "P20_ParentB"="새로미_Parent2"
)


## 2) PCA
pca_fit <- try(prcomp(geno.gwas, center = TRUE, scale. = FALSE), silent = TRUE)
if (inherits(pca_fit, "try-error")) {
  stop("prcomp 실패: geno.gwas가 비었거나 NA/Inf가 과다합니다. 전 단계(Imputation/필터링)를 점검하세요.")
}
scores <- pca_fit$x

## 3) 샘플명 처리 (X 제거 + P매핑)
samples_orig <- rownames(geno.gwas)
if (is.null(samples_orig)) stop("샘플명이 없습니다: rownames(geno.gwas)를 확인하세요.")

samples_show <- sub("^X", "", samples_orig)

is_red <- grepl("^P", samples_show)
samples_mapped <- ifelse(samples_show %in% names(map_P_to_name),
                         unname(map_P_to_name[samples_show]),
                         samples_show)

## 4) 분산 설명력
expvar <- 100 * (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
PCA1 <- round(expvar[1], 3)
PCA2 <- round(expvar[2], 3)

## 5) ggplot용 데이터
df <- data.frame(
  PC1 = scores[,1],
  PC2 = scores[,2],
  sample_show = samples_show,
  label = samples_mapped,
  is_red = is_red,
  color = ifelse(is_red, "red", "gray"),
  stringsAsFactors = FALSE
)

## 6) ggrepel로 겹침 최소화
p <- ggplot(df, aes(PC1, PC2)) +
  geom_point(aes(color = color), size = 2, alpha = 0.85) +
  scale_color_identity() +
  geom_text_repel(
    data = subset(df, is_red),
    aes(label = label),
    color = "red",
    size = 3.2,
    # ✅ 겹침 제어 핵심 파라미터들
    box.padding = 0.35,      # 라벨 상자 간 여백
    point.padding = 0.25,    # 점과 라벨 사이 여백
    force = 1.5,             # 밀어내는 힘(↑면 더 멀리)
    max.overlaps = Inf,      # 라벨이 많아도 포기하지 않음
    min.segment.length = 0,  # 연결선 항상 표시
    segment.color = "red",
    segment.alpha = 0.6,
    seed = 1
  ) +
  labs(
    x = paste0("PC1: ", PCA1, "%"),
    y = paste0("PC2: ", PCA2, "%"),
    title = "PCA Plot (red = samples originally starting with 'P'; labels = mapped names)"
  ) +
  theme_classic(base_family = ifelse(is.na(font_path), "", "korean")) +
  theme(plot.margin = margin(12, 20, 12, 12))

## 7) 저장 (PDF)
out_pdf <- "Figure3_PCA_colored_labelMapped_KR_repel.pdf"
CairoPDF(out_pdf, width = 12.5, height = 8.5)
print(p)
dev.off()

if (!is.na(font_path)) showtext_auto(FALSE)

message("완료: ", out_pdf)


====================================== 
======Windows에서 한글 안깨지고 텍스트 선택 됨========= 
## 1. 필요한 라이브러리 로드 (showtext, Cairo 불필요)
library(dendextend)
library(factoextra)

## ---- 데이터 준비 (기존과 동일) ----
Eucl <- dist(geno.gwas)
fit  <- hclust(Eucl, method = "ward.D2")
dend <- as.dendrogram(fit)

labs_orig <- labels(dend)
labs_show <- sub("^X", "", labs_orig)

# 매핑 정보
map_P_to_name <- c(
  "P01"="굿모닝",
  "P01_VirtualF1"="굿모닝_VirtualF1",
  "P01_ParentA"="굿모닝_Parent1",
  "P01_ParentB"="굿모닝_Parent2",
  "P02"="단청",
  "P02_VirtualF1"="단청_VirtualF1",
  "P02_ParentA"="단청_Parent1",
  "P02_ParentB"="단청_Parent2",
  "P03"="설향",
  "P03_VirtualF1"="설향_VirtualF1",
  "P03_ParentA"="설향_Parent1",
  "P03_ParentB"="설향_Parent2",
  "P04"="봄봄",
  "P04_VirtualF1"="봄봄_VirtualF1",
  "P04_ParentA"="봄봄_Parent1",
  "P04_ParentB"="봄봄_Parent2",  
  "P05"="청춘",
  "P05_VirtualF1"="청춘_VirtualF1",
  "P05_ParentA"="청춘_Parent1",
  "P05_ParentB"="청춘_Parent2",
  "P06"="청어람",
  "P06_VirtualF1"="청어람_VirtualF1",
  "P06_ParentA"="청어람_Parent1",
  "P06_ParentB"="청어람_Parent2",
  "P07"="백리향",
  "P07_VirtualF1"="백리향_VirtualF1",
  "P07_ParentA"="백리향_Parent1",
  "P07_ParentB"="백리향_Parent2",
  "P08"="천리향",
  "P08_VirtualF1"="천리향_VirtualF1",
  "P08_ParentA"="천리향_Parent1",
  "P08_ParentB"="천리향_Parent2",
  "P09"="기라성",
  "P09_VirtualF1"="기라성_VirtualF1",
  "P09_ParentA"="기라성_Parent1",
  "P09_ParentB"="기라성_Parent2",
  "P10"="스트롱",
  "P10_VirtualF1"="스트롱_VirtualF1",
  "P10_ParentA"="스트롱_Parent1",
  "P10_ParentB"="스트롱_Parent2",
  "P11"="베테랑",
  "P11_VirtualF1"="베테랑_VirtualF1",
  "P11_ParentA"="베테랑_Parent1",
  "P11_ParentB"="베테랑_Parent2",
  "P12"="하이탑",
  "P12_VirtualF1"="하이탑_VirtualF1",
  "P12_ParentA"="하이탑_Parent1",
  "P12_ParentB"="하이탑_Parent2",
  "P13"="썸머킹",
  "P13_VirtualF1"="썸머킹_VirtualF1",
  "P13_ParentA"="썸머킹_Parent1",
  "P13_ParentB"="썸머킹_Parent2",
  "P14"="진한백다다기",
  "P14_VirtualF1"="진한백다다기_VirtualF1",
  "P14_ParentA"="진한백다다기_Parent1",
  "P14_ParentB"="진한백다다기_Parent2",
  "P15"="낙동청장",
  "P15_VirtualF1"="낙동청장_VirtualF1",
  "P15_ParentA"="낙동청장_Parent1",
  "P15_ParentB"="낙동청장_Parent2",
  "P16"="신동청장",
  "P16_VirtualF1"="신동청장_VirtualF1",
  "P16_ParentA"="신동청장_Parent1",
  "P16_ParentB"="신동청장_Parent2",
  "P17"="강풍삼척",
  "P17_VirtualF1"="강풍삼척_VirtualF1",
  "P17_ParentA"="강풍삼척_Parent1",
  "P17_ParentB"="강풍삼척_Parent2",
  "P18"="장형낙합",
  "P18_VirtualF1"="장형낙합_VirtualF1",
  "P18_ParentA"="장형낙합_Parent1",
  "P18_ParentB"="장형낙합_Parent2",
  "P19"="신세대",
  "P19_VirtualF1"="신세대_VirtualF1",
  "P19_ParentA"="신세대_Parent1",
  "P19_ParentB"="신세대_Parent2",
  "P20"="새로미",
  "P20_VirtualF1"="새로미_VirtualF1",
  "P20_ParentA"="새로미_Parent1",
  "P20_ParentB"="새로미_Parent2"
)

labs_mapped <- ifelse(
  labs_show %in% names(map_P_to_name),
  unname(map_P_to_name[labs_show]),
  labs_show
)

dend <- set(dend, "labels", labs_mapped)

is_red <- grepl("^P", labs_show)
dend <- dend %>%
 set("labels_col", ifelse(is_red, "red", "black")) %>%
 set("labels_cex", 0.4)

k <- 4
dend_colored <- color_branches(dend, k = k)

## ---- PDF 저장 (Windows 맞춤 설정) ----

# R 내장 함수인 cairo_pdf를 사용합니다.
# family에 윈도우에 무조건 있는 "Malgun Gothic"을 지정합니다.
cairo_pdf("Figure4_dendrogram_Windows2.pdf", 
          width = 13, height = 13, 
          family = "Malgun Gothic")

fviz_dend(
  dend_colored,
  type       = "circular",
  cex        = 0.8,         # 글자 크기가 너무 작으면 0.8 정도로 키워보세요
  label_cols = NULL,
  main       = "Circular Dendrogram (Windows)"
)

dev.off()


=================================== 
## --- 전제: geno.gwas (행=샘플, 열=마커)가 존재 ---
stopifnot(is.matrix(geno.gwas) || is.data.frame(geno.gwas))

## 1) PCA 계산 (center=TRUE, scale=FALSE)
pca_fit <- try(prcomp(geno.gwas, center = TRUE, scale. = FALSE), silent = TRUE)
if (inherits(pca_fit, "try-error")) {
  stop("prcomp 실패: geno.gwas가 비었거나 NA/Inf가 과다합니다. 전 단계(Imputation/필터링)를 점검하세요.")
}

scores <- pca_fit$x  # 샘플 x PC
stopifnot(nrow(scores) == nrow(geno.gwas))

## 2) 샘플명/색상/라벨 선택
samples <- rownames(geno.gwas)
if (is.null(samples)) {
  stop("샘플명이 없습니다: rownames(geno.gwas)를 확인하세요.")
}
is_red  <- grepl("^P", samples)
colors  <- ifelse(is_red, "red", "gray")

## 3) PC1/PC2 설명분산(%) 계산
expvar <- 100 * (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
PCA1 <- round(expvar[1], 3)
PCA2 <- round(expvar[2], 3)

## 4) PDF 장치 열기 & 그리기 (빈 파일 방지용 on.exit)
pdf("Figure3_merged_colored_labeled.pdf", width = 8, height = 6)
on.exit(dev.off(), add = TRUE)

plot(scores[, 1], scores[, 2],
     col = colors, pch = 20, cex = 0.8,
     xlab = paste0("PC1: ", PCA1, "%"),
     ylab = paste0("PC2: ", PCA2, "%"),
     main = "PCA Plot (red_dot=sample name P start; label to red dot)")

## 라벨: 빨간 점만
if (any(is_red)) {
  text(scores[is_red, 1], scores[is_red, 2],
       labels = samples[is_red],
       pos = 4, col = "red", cex = 0.7)
}

legend("topright", legend = c("P_start_samples", "other_samples"),
       col = c("red", "gray"), pch = 20, bty = "n")

dev.off()

========3D PCA================ 
## --- 사용자 제공 매핑 --------------------------------------------------
map_P_to_name <- c(
  "P01"="굿모닝",
  "P01_VirtualF1"="굿모닝_VirtualF1",
  "P01_ParentA"="굿모닝_Parent1",
  "P01_ParentB"="굿모닝_Parent2",
  "P02"="단청",
  "P02_VirtualF1"="단청_VirtualF1",
  "P02_ParentA"="단청_Parent1",
  "P02_ParentB"="단청_Parent2",
  "P03"="설향",
  "P03_VirtualF1"="설향_VirtualF1",
  "P03_ParentA"="설향_Parent1",
  "P03_ParentB"="설향_Parent2",
  "P04"="봄봄",
  "P04_VirtualF1"="봄봄_VirtualF1",
  "P04_ParentA"="봄봄_Parent1",
  "P04_ParentB"="봄봄_Parent2",  
  "P05"="청춘",
  "P05_VirtualF1"="청춘_VirtualF1",
  "P05_ParentA"="청춘_Parent1",
  "P05_ParentB"="청춘_Parent2",
  "P06"="청어람",
  "P06_VirtualF1"="청어람_VirtualF1",
  "P06_ParentA"="청어람_Parent1",
  "P06_ParentB"="청어람_Parent2",
  "P07"="백리향",
  "P07_VirtualF1"="백리향_VirtualF1",
  "P07_ParentA"="백리향_Parent1",
  "P07_ParentB"="백리향_Parent2",
  "P08"="천리향",
  "P08_VirtualF1"="천리향_VirtualF1",
  "P08_ParentA"="천리향_Parent1",
  "P08_ParentB"="천리향_Parent2",
  "P09"="기라성",
  "P09_VirtualF1"="기라성_VirtualF1",
  "P09_ParentA"="기라성_Parent1",
  "P09_ParentB"="기라성_Parent2",
  "P10"="스트롱",
  "P10_VirtualF1"="스트롱_VirtualF1",
  "P10_ParentA"="스트롱_Parent1",
  "P10_ParentB"="스트롱_Parent2",
  "P11"="베테랑",
  "P11_VirtualF1"="베테랑_VirtualF1",
  "P11_ParentA"="베테랑_Parent1",
  "P11_ParentB"="베테랑_Parent2",
  "P12"="하이탑",
  "P12_VirtualF1"="하이탑_VirtualF1",
  "P12_ParentA"="하이탑_Parent1",
  "P12_ParentB"="하이탑_Parent2",
  "P13"="썸머킹",
  "P13_VirtualF1"="썸머킹_VirtualF1",
  "P13_ParentA"="썸머킹_Parent1",
  "P13_ParentB"="썸머킹_Parent2",
  "P14"="진한백다다기",
  "P14_VirtualF1"="진한백다다기_VirtualF1",
  "P14_ParentA"="진한백다다기_Parent1",
  "P14_ParentB"="진한백다다기_Parent2",
  "P15"="낙동청장",
  "P15_VirtualF1"="낙동청장_VirtualF1",
  "P15_ParentA"="낙동청장_Parent1",
  "P15_ParentB"="낙동청장_Parent2",
  "P16"="신동청장",
  "P16_VirtualF1"="신동청장_VirtualF1",
  "P16_ParentA"="신동청장_Parent1",
  "P16_ParentB"="신동청장_Parent2",
  "P17"="강풍삼척",
  "P17_VirtualF1"="강풍삼척_VirtualF1",
  "P17_ParentA"="강풍삼척_Parent1",
  "P17_ParentB"="강풍삼척_Parent2",
  "P18"="장형낙합",
  "P18_VirtualF1"="장형낙합_VirtualF1",
  "P18_ParentA"="장형낙합_Parent1",
  "P18_ParentB"="장형낙합_Parent2",
  "P19"="신세대",
  "P19_VirtualF1"="신세대_VirtualF1",
  "P19_ParentA"="신세대_Parent1",
  "P19_ParentB"="신세대_Parent2",
  "P20"="새로미",
  "P20_VirtualF1"="새로미_VirtualF1",
  "P20_ParentA"="새로미_Parent1",
  "P20_ParentB"="새로미_Parent2"
)

## --- 샘플ID 정리(공백 제거) --------------------------------------------
scores$Sample <- trimws(as.character(scores$Sample))

## --- 실제 존재하는 샘플만 하이라이트 대상으로 --------------------------
highlight_samples <- intersect(names(map_P_to_name), scores$Sample)

cat("하이라이트 대상(실제 데이터에 존재):", length(highlight_samples), "\n")
if (length(setdiff(names(map_P_to_name), scores$Sample)) > 0) {
  cat("매핑에는 있지만 데이터에는 없는 샘플(무시됨):\n")
  print(setdiff(names(map_P_to_name), scores$Sample))
}

## --- 라벨/색/크기 컬럼 생성 --------------------------------------------
scores$Highlight <- scores$Sample %in% highlight_samples

scores$Label <- ifelse(scores$Highlight,
                       unname(map_P_to_name[scores$Sample]),
                       scores$Sample)

scores$Color <- ifelse(scores$Highlight, "red", "grey70")
scores$Size  <- ifelse(scores$Highlight, 6, 3)

png("PCA_2D_PC1_PC2_highlightP.png", width = 1000, height = 900, res = 120)

plot(scores$PC1, scores$PC2,
     pch = 19, col = scores$Color, cex = scores$Size/3,
     xlab = paste0("PC1 (", pc1p, "%)"),
     ylab = paste0("PC2 (", pc2p, "%)"),
     main = "PCA (PC1 vs PC2) - Highlight P samples")

idx <- which(scores$Highlight)
text(scores$PC1[idx], scores$PC2[idx],
     labels = scores$Label[idx],
     col = "red4", cex = 0.8, pos = 3)

legend("topright",
       legend = c("Highlighted(P-set)", "Others"),
       col = c("red", "grey70"), pch = 19, pt.cex = c(1.4, 1),
       bty = "n")

dev.off()


suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmlwidgets))

p3d <- plot_ly(
  data = scores,
  x = ~PC1, y = ~PC2, z = ~PC3,
  type = "scatter3d", mode = "markers",
  marker = list(color = ~Color, size = ~Size, opacity = 0.85),
  text = ~paste0(Label, " [", Sample, "]"),
  hoverinfo = "text"
) %>%
  layout(
    title = "PCA 3D (PC1-PC2-PC3) - Highlight P samples",
    scene = list(
      xaxis = list(title = paste0("PC1 (", pc1p, "%)")),
      yaxis = list(title = paste0("PC2 (", pc2p, "%)")),
      zaxis = list(title = paste0("PC3 (", pc3p, "%)"))
    )
  )

p3d
saveWidget(p3d, "PCA_3D_PC1_PC2_PC3_highlightP.html", selfcontained = TRUE)

======================= 
###plotly 3D에서 “텍스트가 항상 보이게” 만들기 (핵심)
library(plotly)
library(htmlwidgets)

# 하이라이트/라벨 데이터 분리
df_hi  <- scores[scores$Highlight, ]
df_oth <- scores[!scores$Highlight, ]

p3d <- plot_ly() %>%
  # 1) 나머지 샘플: 점만 + hover로 라벨
  add_trace(
    data = df_oth,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "grey70", opacity = 0.65),
    text = ~paste0(Sample),
    hoverinfo = "text",
    name = "Others"
  ) %>%
  # 2) 하이라이트 샘플: 점 + 텍스트(항상 보이게)
  add_trace(
    data = df_hi,
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers+text",
    marker = list(size = 6, color = "red", opacity = 0.95),
    text = ~paste0(Label),              # ← “주어진 샘플명(품종명)” 표시
    textposition = "top center",
    textfont = list(color = "red4", size = 12),
    hoverinfo = "text",
    hovertext = ~paste0(Label, " [", Sample, "]"),
    name = "Highlighted"
  ) %>%
  layout(
    title = "PCA 3D (PC1-PC2-PC3) - Highlight labels",
    scene = list(
      xaxis = list(title = paste0("PC1 (", pc1p, "%)")),
      yaxis = list(title = paste0("PC2 (", pc2p, "%)")),
      zaxis = list(title = paste0("PC3 (", pc3p, "%)"))
    ),
    legend = list(orientation = "h", x = 0.05, y = 1.02)
  )

p3d
saveWidget(p3d, "PCA_3D_highlight_with_labels.html", selfcontained = TRUE)
#file:///C:/Users/2230075/AppData/Local/Temp/RtmpyM4Ak5/viewhtml459c6737ff8/index.html



================================ 
## (A) scores_df 생성 라인 오타 수정
# 기존 잘못된 코드:
# scores_df <- .data.frame(scores, stringsAsFactors = FALSE)

# 올바른 코드:
scores_df <- as.data.frame(scores, stringsAsFactors = FALSE)

# prcomp 이름 보정(혹시를 대비해)
colnames(scores_df) <- paste0("PC", seq_len(ncol(scores_df)))
scores_df$Sample <- samples

## (B) 나머지 파이프라인 계속 진행
expvar <- 100 * (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
pc1p <- if (length(expvar) >= 1) round(expvar[1], 3) else NA_real_
pc2p <- if (length(expvar) >= 2) round(expvar[2], 3) else NA_real_
pc3p <- if (length(expvar) >= 3) round(expvar[3], 3) else NA_real_

map_P_to_name <- c(
  "P01"="굿모닝","P01_VirtualF1"="굿모닝_VirtualF1","P01_ParentA"="굿모닝_Parent1","P01_ParentB"="굿모닝_Parent2",
  "P02"="단청","P02_VirtualF1"="단청_VirtualF1","P02_ParentA"="단청_Parent1","P02_ParentB"="단청_Parent2",
  "P03"="설향","P03_VirtualF1"="설향_VirtualF1","P03_ParentA"="설향_Parent1","P03_ParentB"="설향_Parent2",
  "P04"="봄봄","P04_VirtualF1"="봄봄_VirtualF1","P04_ParentA"="봄봄_Parent1","P04_ParentB"="봄봄_Parent2",
  "P05"="청춘","P05_VirtualF1"="청춘_VirtualF1","P05_ParentA"="청춘_Parent1","P05_ParentB"="청춘_Parent2",
  "P06"="청어람","P06_VirtualF1"="청어람_VirtualF1","P06_ParentA"="청어람_Parent1","P06_ParentB"="청어람_Parent2",
  "P07"="백리향","P07_VirtualF1"="백리향_VirtualF1","P07_ParentA"="백리향_Parent1","P07_ParentB"="백리향_Parent2",
  "P08"="천리향","P08_VirtualF1"="천리향_VirtualF1","P08_ParentA"="천리향_Parent1","P08_ParentB"="천리향_Parent2",
  "P09"="기라성","P09_VirtualF1"="기라성_VirtualF1","P09_ParentA"="기라성_Parent1","P09_ParentB"="기라성_Parent2",
  "P10"="스트롱","P10_VirtualF1"="스트롱_VirtualF1","P10_ParentA"="스트롱_Parent1","P10_ParentB"="스트롱_Parent2",
  "P11"="베테랑","P11_VirtualF1"="베테랑_VirtualF1","P11_ParentA"="베테랑_Parent1","P11_ParentB"="베테랑_Parent2",
  "P12"="하이탑","P12_VirtualF1"="하이탑_VirtualF1","P12_ParentA"="하이탑_Parent1","P12_ParentB"="하이탑_Parent2",
  "P13"="썸머킹","P13_VirtualF1"="썸머킹_VirtualF1","P13_ParentA"="썸머킹_Parent1","P13_ParentB"="썸머킹_Parent2",
  "P14"="진한백다다기","P14_VirtualF1"="진한백다다기_VirtualF1","P14_ParentA"="진한백다다기_Parent1","P14_ParentB"="진한백다다기_Parent2",
  "P15"="낙동청장","P15_VirtualF1"="낙동청장_VirtualF1","P15_ParentA"="낙동청장_Parent1","P15_ParentB"="낙동청장_Parent2",
  "P16"="신동청장","P16_VirtualF1"="신동청장_VirtualF1","P16_ParentA"="신동청장_Parent1","P16_ParentB"="신동청장_Parent2",
  "P17"="강풍삼척","P17_VirtualF1"="강풍삼척_VirtualF1","P17_ParentA"="강풍삼척_Parent1","P17_ParentB"="강풍삼척_Parent2",
  "P18"="장형낙합","P18_VirtualF1"="장형낙합_VirtualF1","P18_ParentA"="장형낙합_Parent1","P18_ParentB"="장형낙합_Parent2",
  "P19"="신세대","P19_VirtualF1"="신세대_VirtualF1","P19_ParentA"="신세대_Parent1","P19_ParentB"="신세대_Parent2",
  "P20"="새로미","P20_VirtualF1"="새로미_VirtualF1","P20_ParentA"="새로미_Parent1","P20_ParentB"="새로미_Parent2"
)

highlight_samples <- intersect(names(map_P_to_name), scores_df$Sample)

scores_df$Highlight <- scores_df$Sample %in% highlight_samples
scores_df$Label     <- ifelse(scores_df$Highlight,
                              unname(map_P_to_name[scores_df$Sample]),
                              scores_df$Sample)
scores_df$Color     <- ifelse(scores_df$Highlight, "red", "grey70")
scores_df$Size      <- ifelse(scores_df$Highlight, 6, 3)

## 2D (PC1 vs PC2)
png("PCA_2D_PC1_PC2_highlightP.png", width = 1000, height = 900, res = 120)
plot(scores_df$PC1, scores_df$PC2,
     pch = 19, col = scores_df$Color, cex = scores_df$Size/3,
     xlab = paste0("PC1 (", pc1p, "%)"),
     ylab = paste0("PC2 (", pc2p, "%)"),
     main = "PCA (PC1 vs PC2) - Highlight P samples")
idx <- which(scores_df$Highlight)
if (length(idx) > 0) {
  text(scores_df$PC1[idx], scores_df$PC2[idx],
       labels = scores_df$Label[idx],
       col = "red4", cex = 0.8, pos = 3)
}
legend("topright",
       legend = c("Highlighted(P-set)", "Others"),
       col = c("red", "grey70"), pch = 19, pt.cex = c(1.4, 1),
       bty = "n")
dev.off()

## 3D (PC1-PC2-PC3) — PC3가 있으면만
if (ncol(scores_df) >= 4) {
  suppressPackageStartupMessages(library(plotly))
  suppressPackageStartupMessages(library(htmlwidgets))
  df_hi  <- subset(scores_df, Highlight)
  df_oth <- subset(scores_df, !Highlight)

  p3d <- plot_ly() %>%
    add_trace(data = df_oth, x = ~PC1, y = ~PC2, z = ~PC3,
              type = "scatter3d", mode = "markers",
              marker = list(size = 3, color = "grey70", opacity = 0.65),
              text = ~Sample, hoverinfo = "text", name = "Others") %>%
    add_trace(data = df_hi, x = ~PC1, y = ~PC2, z = ~PC3,
              type = "scatter3d", mode = "markers+text",
              marker = list(size = 6, color = "red", opacity = 0.95),
              text = ~Label, textposition = "top center",
              textfont = list(color = "red4", size = 12),
              hoverinfo = "text",
              hovertext = ~paste0(Label, " [", Sample, "]"),
              name = "Highlighted") %>%
    layout(
      title = "PCA 3D (PC1-PC2-PC3) - Highlight labels",
      scene = list(
        xaxis = list(title = paste0("PC1 (", pc1p, "%)")),
        yaxis = list(title = paste0("PC2 (", pc2p, "%)")),
        zaxis = list(title = paste0("PC3 (", pc3p, "%)"))
      ),
      legend = list(orientation = "h", x = 0.05, y = 1.02)
    )

  saveWidget(p3d, "PCA_3D_highlight_with_labels.html", selfcontained = TRUE)
} else {
  message("PC3가 없어 3D plotly는 생략되었습니다. (주성분 수: ", ncol(scores_df) - 1, ")")
}



=========================== 
## 0) 패키지
suppressPackageStartupMessages({
  library(dendextend)
  library(factoextra)
  library(ggplot2)
  library(ggforce)    # 군집 바운더리(말풍선형/볼록헐)용
})

## 1) 거리/군집/덴드로그램
Eucl <- dist(geno.gwas)
fit  <- hclust(Eucl, method = "ward.D2")
dend <- as.dendrogram(fit)

## 2) 라벨 매핑(사용자 제공)
labs_orig <- labels(dend)
labs_show <- sub("^X", "", labs_orig)

map_P_to_name <- c(
  "P01"="굿모닝","P01_VirtualF1"="굿모닝_VirtualF1","P01_ParentA"="굿모닝_Parent1","P01_ParentB"="굿모닝_Parent2",
  "P02"="단청","P02_VirtualF1"="단청_VirtualF1","P02_ParentA"="단청_Parent1","P02_ParentB"="단청_Parent2",
  "P03"="설향","P03_VirtualF1"="설향_VirtualF1","P03_ParentA"="설향_Parent1","P03_ParentB"="설향_Parent2",
  "P04"="봄봄","P04_VirtualF1"="봄봄_VirtualF1","P04_ParentA"="봄봄_Parent1","P04_ParentB"="봄봄_Parent2",
  "P05"="청춘","P05_VirtualF1"="청춘_VirtualF1","P05_ParentA"="청춘_Parent1","P05_ParentB"="청춘_Parent2",
  "P06"="청어람","P06_VirtualF1"="청어람_VirtualF1","P06_ParentA"="청어람_Parent1","P06_ParentB"="청어람_Parent2",
  "P07"="백리향","P07_VirtualF1"="백리향_VirtualF1","P07_ParentA"="백리향_Parent1","P07_ParentB"="백리향_Parent2",
  "P08"="천리향","P08_VirtualF1"="천리향_VirtualF1","P08_ParentA"="천리향_Parent1","P08_ParentB"="천리향_Parent2",
  "P09"="기라성","P09_VirtualF1"="기라성_VirtualF1","P09_ParentA"="기라성_Parent1","P09_ParentB"="기라성_Parent2",
  "P10"="스트롱","P10_VirtualF1"="스트롱_VirtualF1","P10_ParentA"="스트롱_Parent1","P10_ParentB"="스트롱_Parent2",
  "P11"="베테랑","P11_VirtualF1"="베테랑_VirtualF1","P11_ParentA"="베테랑_Parent1","P11_ParentB"="베테랑_Parent2",
  "P12"="하이탑","P12_VirtualF1"="하이탑_VirtualF1","P12_ParentA"="하이탑_Parent1","P12_ParentB"="하이탑_Parent2",
  "P13"="썸머킹","P13_VirtualF1"="썸머킹_VirtualF1","P13_ParentA"="썸머킹_Parent1","P13_ParentB"="썸머킹_Parent2",
  "P14"="진한백다다기","P14_VirtualF1"="진한백다다기_VirtualF1","P14_ParentA"="진한백다다기_Parent1","P14_ParentB"="진한백다다기_Parent2",
  "P15"="낙동청장","P15_VirtualF1"="낙동청장_VirtualF1","P15_ParentA"="낙동청장_Parent1","P15_ParentB"="낙동청장_Parent2",
  "P16"="신동청장","P16_VirtualF1"="신동청장_VirtualF1","P16_ParentA"="신동청장_Parent1","P16_ParentB"="신동청장_Parent2",
  "P17"="강풍삼척","P17_VirtualF1"="강풍삼척_VirtualF1","P17_ParentA"="강풍삼척_Parent1","P17_ParentB"="강풍삼척_Parent2",
  "P18"="장형낙합","P18_VirtualF1"="장형낙합_VirtualF1","P18_ParentA"="장형낙합_Parent1","P18_ParentB"="장형낙합_Parent2",
  "P19"="신세대","P19_VirtualF1"="신세대_VirtualF1","P19_ParentA"="신세대_Parent1","P19_ParentB"="신세대_Parent2",
  "P20"="새로미","P20_VirtualF1"="새로미_VirtualF1","P20_ParentA"="새로미_Parent1","P20_ParentB"="새로미_Parent2"
)

labs_mapped <- ifelse(labs_show %in% names(map_P_to_name),
                      unname(map_P_to_name[labs_show]),
                      labs_show)

## 덴드에 표시 라벨 적용 및 기본 스타일
is_red <- grepl("^P", labs_show)
dend <- dend %>%
  set("labels", labs_mapped) %>%
  set("labels_col", ifelse(is_red, "red", "black")) %>%
  set("labels_cex", 0.4)

## 3) 군집(4개) 색칠
k <- 4
dend_colored <- color_branches(dend, k = k)

## 4) 원형 덴드로그램 ggplot 객체
p <- fviz_dend(
  dend_colored,
  type       = "circular",
  cex        = 0.8,
  label_cols = NULL,
  main       = "Circular Dendrogram with Group Hulls"
)

## 5) 리프 좌표 추출 → 군집/그룹 매핑
# cutree는 '원래 샘플ID(=labs_show)' 기준으로 군집 번호를 반환합니다.
cl_vec <- cutree(fit, k = k)   # named integer vector (names = 원샘플)
# fviz_dend가 만든 ggplot에서 라벨 레이어 데이터만 꺼냅니다.
gb  <- ggplot_build(p)
lab_layers <- gb$data[ vapply(gb$data, function(d) "label" %in% names(d), logical(1)) ]
stopifnot(length(lab_layers) >= 1)
lbl <- lab_layers[[length(lab_layers)]]   # 마지막 레이어에 라벨이 있는 경우가 보통

# 표시라벨(=labs_mapped) -> 원샘플ID(=labs_show) 역매핑
display_to_id <- setNames(labs_show, labs_mapped)
lbl$SampleID   <- display_to_id[as.character(lbl$label)]
lbl$Cluster    <- cl_vec[ lbl$SampleID ]

## 6) 군집번호 -> 그룹명(원하시는 이름으로 지정)
# ⚠️ cutree 번호는 데이터마다 달라질 수 있습니다.
# 아래 출력으로 각 군집에 어떤 샘플이 들어있는지 확인 후, 이름을 맞춰주세요.
cat("\n=== 군집별 예시 라벨(앞 5개) ===\n")
for (i in sort(unique(lbl$Cluster))) {
  cat(sprintf("Cluster %s : %s\n", i,
      paste(head(lbl$label[lbl$Cluster==i], 5), collapse = ", ")))
}

# 기본 매핑(필요시 숫자를 바꾸세요)
cluster_to_group <- c(
  `1` = "백다다기",  # khaki
  `2` = "취청오이",  # navy
  `3` = "가시오이",  # purple
  `4` = "미니오이"   # red
)

lbl$Group <- factor(cluster_to_group[as.character(lbl$Cluster)],
                    levels = c("백다다기","취청오이","가시오이","미니오이"))

## 7) 그룹 색상(바운더리 색/채우기)
group_cols <- c(
  "백다다기" = "#B9A46C",  # khaki-like
  "취청오이" = "#19398A",  # navy
  "가시오이" = "#6A0DAD",  # purple
  "미니오이" = "#D62728"   # red
)

## 8) 바운더리(볼록헐/말풍선) + 그룹명 라벨 추가
# - geom_mark_hull: 그룹을 감싸는 폴리곤과 그룹명 레이블을 동시에 표현
# - expand/label.* 파라미터로 여백과 라벨 스타일 조절
p_hull <- p +
  ggforce::geom_mark_hull(
    data = lbl,
    aes(x, y, group = Group, fill = Group, label = Group),
    alpha = 0.1,        # 투명도
    con.cap = unit(0, "mm"),
    con.size = 0.2,
    expand = unit(6, "mm"),
    label.fontsize = 12,
    label.buffer = unit(3, "mm"),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_cols)

## 9) 저장(Windows: Malgun Gothic 내장 장점)
ggsave("Figure4_dendrogram_with_group_hulls.pdf",
       plot = p_hull, device = cairo_pdf, width = 13, height = 13, family = "Malgun Gothic")
ggsave("Figure4_dendrogram_with_group_hulls.png",
       plot = p_hull, width = 2400, height = 2400, units = "px", dpi = 300)

================================== 
# 그룹 중심을 계산해 라벨을 바깥쪽으로 이동(이미 앞선 코드에 유사 로직 포함)
group_centers <- aggregate(cbind(x, y) ~ Group, data = lbl, FUN = mean)
group_centers$r  <- sqrt(group_centers$x^2 + group_centers$y^2)
group_centers$xf <- group_centers$x / group_centers$r * group_centers$r * 1.10
group_centers$yf <- group_centers$y / group_centers$r * group_centers$r * 1.10

p_group <- p_under +
  geom_text(
    data = group_centers,
    aes(x = xf, y = yf, label = Group),
    colour = "grey20", fontface = "bold",
    size = 0.0, show.legend = FALSE
  )

ggsave("Figure4_dendrogram_group_only.pdf",
       plot = p_group, device = cairo_pdf, width = 13, height = 13,
       family = "Malgun Gothic")

======================== 
library(ggplot2)

# 0) 말굽고딕(맑은 고딕) 사용 중이면 이후 ggsave에서 family="Malgun Gothic"와 매칭됩니다.
# windowsFonts(MalgunGothic = windowsFont("Malgun Gothic"))

# 1) p_under 복구 시도
if (exists("p_under") && inherits(p_under, "ggplot")) {
  message("p_under 이미 존재합니다.")
} else {
  # last_plot()로 직전 플롯 복구 시도
  tmp_last <- try(ggplot2::last_plot(), silent = TRUE)
  if (!inherits(tmp_last, "try-error") && inherits(tmp_last, "ggplot")) {
    p_under <- tmp_last
    message("last_plot()에서 p_under 복구했습니다.")
  } else {
    message("p_under가 없어 기본 플롯을 새로 생성합니다.")
    # ---- [대체용 최소 플롯] ----
    # lbl: x, y, Group, label을 가진 데이터프레임이라고 가정
    # 원형 배치 좌표(lbl$x, lbl$y)를 이미 갖고 있다면 점과 텍스트만 얹습니다.
    p_under <- ggplot(lbl, aes(x = x, y = y)) +
      geom_point(aes(color = Group), size = 1, alpha = 0.8, show.legend = FALSE) +
      # 원래 원형 dendrogram의 서클/링/가이드가 있다면 여기서 추가하세요.
      coord_equal() +
      theme_void() +
      theme(
        text = element_text(family = "Malgun Gothic"),
        plot.margin = margin(10, 10, 10, 10)
      )
  }
}

# 그룹 중심 계산
group_centers <- aggregate(cbind(x, y) ~ Group, data = lbl, FUN = mean)
group_centers$r  <- sqrt(group_centers$x^2 + group_centers$y^2)

# 바깥으로 밀어낼 배율 (필요시 1.10 -> 1.20~1.30로 키우세요)
push_factor <- 1.20
group_centers$xf <- group_centers$x / group_centers$r * group_centers$r * push_factor
group_centers$yf <- group_centers$y / group_centers$r * group_centers$r * push_factor

# 그룹 라벨 추가 (텍스트가 위로 올라오도록 마지막 레이어에 추가)
p_group <- p_under +
  geom_text(
    data = group_centers,
    aes(x = xf, y = yf, label = Group),
    colour = "grey20", fontface = "bold",
    size = 4,            # ← 0이 아니어야 보입니다. (필요시 조절)
    show.legend = FALSE
  )
  
  
ggsave(
  filename = "Figure4_dendrogram_group_only.pdf",
  plot = p_group,
  device = cairo_pdf,
  width = 13, height = 13,
  family = "Malgun Gothic"  # 윈도우에서 맑은 고딕
)

======================= 
## ----  A. 필요 패키지 로드 (추가) ----
# install.packages(c("ape", "dplyr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ggtree")

library(ape)
library(dplyr)
library(ggtree)

## ----  B. hclust -> phylo 변환 & tip 라벨을 'labs_mapped'로 통일 ----
phy <- as.phylo(fit)

# fit/dend와 동일한 라벨 매핑 테이블: labs_orig -> labs_mapped
# (labs_orig/labs_mapped은 사용자 코드에서 이미 계산됨)
map_label <- setNames(labs_mapped, labs_orig)

# phy의 tip.label을 사용자가 만든 매핑으로 치환
tmp <- unname(map_label[phy$tip.label])
# 매핑이 없는 항목(KXXXX 등)은 X 접두어만 제거하고 그대로 사용
phy$tip.label <- ifelse(is.na(tmp), sub("^X", "", phy$tip.label), tmp)

## ----  C. 4개 그룹으로 분할 (cutree) & 라벨을 phy와 동일 체계로 정렬 ----
# cutree 결과는 원래 행이름(= labs_orig) 기준이므로, labs_mapped로 이름을 바꿔 phy와 일치시킴
grp_raw <- cutree(fit, k = k)  # k = 4 (상단 코드에서 지정)
names(grp_raw) <- names(grp_raw)                    # 방어적
# 먼저 X를 제거
names(grp_raw) <- sub("^X", "", names(grp_raw))
# 그 다음, P코드 등 매핑 가능한 것은 한글명으로 변환
names(grp_raw) <- ifelse(
  names(grp_raw) %in% names(map_P_to_name),
  unname(map_P_to_name[names(grp_raw)]),
  names(grp_raw)
)

# 이제 names(grp_raw)와 phy$tip.label이 같은 체계(매핑된 라벨)여야 함
# 그룹 레벨(1~4)을 정렬해 고정된 색상 매핑을 보장
grp_levels <- sort(unique(grp_raw))
tips_by_group <- lapply(grp_levels, function(g) names(grp_raw)[grp_raw == g])
names(tips_by_group) <- grp_levels

# 각 그룹의 MRCA 노드(최소 공통 조상) 계산
mrca_nodes <- sapply(tips_by_group, function(tips) getMRCA(phy, tips))

## ----  D. 색상/스타일(원하는 팔레트로 교체 가능) ----
fill_cols  <- c("#7DB9B6", "#FDAF75", "#A3C4F3", "#C3DFA6")  # 그룹 1~4
border_col <- "grey55"     # 그룹 경계선을 검정보다 연한 그레이로
alpha_fill <- 0.35

# 텍스트/오프셋 설정
label_size   <- 2.6        # ggtree의 포인트 크기 기준
label_offset <- 0.8        # 원형에서 바깥으로 밀기 (겹침 완화)
key_labels   <- c("백다다기","취청오이","가시오이","미니오이")  # 더 멀리/굵게

## ----  E. 원형 덴드로그램 + 곡선 하이라이트(geom_hilight) ----
p_circ <- ggtree(phy, layout = "circular", size = 0.3, color = "grey25")

# 하이라이트(먼저 그린다 → z-order상 아래로 깔림)
for (i in seq_along(mrca_nodes)) {
  p_circ <- p_circ + geom_hilight(
    node   = mrca_nodes[i],
    fill   = fill_cols[i],
    alpha  = alpha_fill,
    colour = border_col
  )
}

# 기본 라벨(마지막에 다시 그려서 텍스트가 위에 올라오게 함)
p_circ <- p_circ + geom_tiplab(
  aes(label = label),
  size   = label_size,
  offset = label_offset,
  color  = "black"
)

# 핵심 라벨 4개만 더 바깥으로, 겹침 최소화를 위해 굵게
lab_pt <- p_circ$data %>% filter(isTip, label %in% key_labels)
p_circ <- p_circ + geom_tiplab(
  data   = lab_pt,
  aes(label = label),
  size   = label_size + 0.4,
  offset = label_offset + 0.8,  # 더 멀리
  fontface = "bold",
  color  = "black"
)

# (선택) 전체 텍스트에 윈도우 폰트 적용
p_circ <- p_circ + theme(
  text = element_text(family = "Malgun Gothic")
)

## ----  F. PDF 저장 (현재 방식 유지: Cairo 내장 글꼴 임베드) ----
ggsave(
  filename = "Figure4_dendrogram_highlight_circ.pdf",
  plot     = p_circ,
  width    = 12, height = 12,
  device   = cairo_pdf
)

## ----  G. 직사각형 레이아웃 버전도 필요하면 참고 ----
# p_rect <- ggtree(phy, layout = "rectangular", size = 0.3, color = "grey25")
# for (i in seq_along(mrca_nodes)) {
#   p_rect <- p_rect + geom_hilight(
#     node   = mrca_nodes[i],
#     fill   = fill_cols[i],
#     alpha  = alpha_fill,
#     colour = border_col
#   )
# }
# p_rect <- p_rect + geom_tiplab(aes(label = label),
#                                size = 3.0, hjust = 0, offset = 0.02, color = "black") +
#   xlim(0, max(p_rect$data$x) * 1.08) +
#   theme(text = element_text(family = "Malgun Gothic"))
# ggsave("Figure4_dendrogram_highlight_rect.pdf", p_rect, width = 13, height = 10, device = cairo_pdf)
  
  
========================= 
## ===================== 0) 패키지 & 폰트 =====================
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2)
  library(ggtree); library(showtext); library(sysfonts)
})

# The most robust method: Download Nanum Gothic directly from Google Fonts.
# This completely bypasses Windows/Mac local font path issues.
tryCatch({
  font_add_google("Nanum Gothic", "nanum")
}, error = function(e) {
  # Fallback just in case you have no internet connection
  if (file.exists("C:/Windows/Fonts/malgun.ttf")) {
    font_add("nanum", "C:/Windows/Fonts/malgun.ttf")
  } else if (file.exists("/System/Library/Fonts/Supplemental/AppleGothic.ttf")) {
    font_add("nanum", "/System/Library/Fonts/Supplemental/AppleGothic.ttf")
  } else {
    message("Warning: No Korean font found!")
  }
})

# Turn on showtext to convert text into shapes (prevents PDF breaking)
showtext_auto(TRUE)
options(ggtree_family = "nanum") # Force ggtree to use the loaded font

## ===================== 1) phy 생성 & 라벨 매핑 =====================
phy <- as.phylo(fit)
map_label <- setNames(labs_mapped, labs_orig)

tmp <- unname(map_label[phy$tip.label])
phy$tip.label <- ifelse(is.na(tmp), sub("^X", "", phy$tip.label), tmp)

## 그룹 만들기
grp_raw <- cutree(fit, k = k)
names(grp_raw) <- sub("^X","", names(grp_raw))
names(grp_raw) <- ifelse(
  names(grp_raw) %in% names(map_P_to_name),
  unname(map_P_to_name[names(grp_raw)]),
  names(grp_raw)
)

grp_levels <- sort(unique(grp_raw))
tips_by_group <- lapply(grp_levels, function(g) names(grp_raw)[grp_raw == g])
names(tips_by_group) <- grp_levels

# MRCA 계산(1개 팁인 경우 대비)
mrca_nodes <- sapply(tips_by_group, function(tips) if (length(tips) < 2) NA_integer_ else getMRCA(phy, tips))

## ===================== 2) 팔레트/스타일 =====================
fill_cols <- c("#7DB9B6", "#FDAF75", "#A3C4F3", "#C3DFA6")
if (length(grp_levels) != length(fill_cols)) {
  fill_cols <- scales::hue_pal()(length(grp_levels))
}
names(fill_cols) <- grp_levels

# 균일 두께 링을 원하는 폭으로:
ring_outer_pad <- 0.6   # 트리 바깥 여백
ring_thick     <- 1.0   # 링 두께 (모든 그룹 동일)

label_size   <- 2.6
label_offset <- 0.9
key_labels   <- c("백다다기","취청오이","가시오이","미니오이")

## ===================== 3) 트리(선)부터: groupOTU + geom_tree =====================
phy_grouped <- groupOTU(phy, tips_by_group)
p_circ <- ggtree(phy_grouped, layout = "circular", aes(color = group)) +
  geom_tree(size = 0.35) +
  scale_color_manual(values = fill_cols, na.value = "transparent") +
  theme_minimal(base_family = getOption("ggtree_family")) +
  theme(legend.position = "none", plot.margin = margin(8,8,8,8))

## ===================== 4) 하이라이트: 버전 호환 처리 =====================
# 현재 ggtree의 geom_hilight가 'offset'을 지원하는지 점검
hl_formals <- names(formals(ggtree::geom_hilight))
support_offset <- "offset" %in% hl_formals
support_extend <- "extend" %in% hl_formals || "extendto" %in% hl_formals

if (support_offset) {
  # --- A안: offset/extend 있는 최신 계열 ---
  hl_offset <- 0.5
  hl_extend <- ring_thick
  for (g in names(mrca_nodes)) {
    nd <- mrca_nodes[[g]]; if (is.na(nd)) next
    p_circ <- p_circ +
      geom_hilight(
        node   = nd,
        fill   = fill_cols[g],
        colour = fill_cols[g],
        alpha  = 0.35,
        offset = hl_offset,
        extend = hl_extend
      )
  }
} else {
  # --- B안: offset 미지원 → 균일 두께 링을 수동 합성 ---
  # 트리 바깥 반지름 계산
  d_tree <- p_circ$data
  R_tip  <- max(d_tree$x, na.rm = TRUE)
  R_out  <- R_tip + ring_outer_pad
  for (g in names(mrca_nodes)) {
    nd <- mrca_nodes[[g]]; if (is.na(nd)) next
    pos <- get_clade_position(p_circ, nd)  # ymin~ymax: 각도 구간, xmin/xmax: 반지름 구간
    p_circ <- p_circ +
      annotate("rect",
               xmin = R_out - ring_thick, xmax = R_out,
               ymin = pos$ymin,          ymax = pos$ymax,
               fill = fill_cols[g], color = fill_cols[g], alpha = 0.35)
  }
  # 라벨을 위한 x축 확장
  p_circ <- p_circ + xlim(0, R_out + label_offset + 0.5)
}

## ===================== 5) 라벨: 검정 → P로 시작하는 것만 빨강 → 핵심 4개 더 바깥 =====================
p_circ <- p_circ +
  geom_tiplab(aes(label = label),
              size = label_size, offset = label_offset,
              color = "black", family = getOption("ggtree_family"))

tip_df  <- p_circ$data %>% filter(isTip)
tip_red <- tip_df %>% filter(grepl("^[pP]", label))   # 필요시 정규식 변경
if (nrow(tip_red) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_red, aes(label = label),
                size = label_size, offset = label_offset,
                color = "red3", fontface = "bold",
                family = getOption("ggtree_family"))
}

lab_pt <- tip_df %>% filter(label %in% key_labels)
if (nrow(lab_pt) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = lab_pt, aes(label = label),
                size = label_size + 0.5,
                offset = label_offset + 0.8,
                fontface = "bold", color = "black",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
ggsave("Figure4_dendrogram_highlight_circ_v3_1.pdf", p_circ,
       width = 12, height = 12, device = cairo_pdf)
ggsave("Figure4_dendrogram_highlight_circ_v3_1.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)
	   
	   
========================================
## ===================== 0) 패키지 & 폰트 =====================
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2)
  library(ggtree); library(showtext); library(sysfonts)
})

# The most robust method for PDFs: Download Nanum Gothic directly from Google Fonts.
tryCatch({
  font_add_google("Nanum Gothic", "nanum")
}, error = function(e) {
  # Fallback if no internet
  if (file.exists("C:/Windows/Fonts/malgun.ttf")) {
    font_add("nanum", "C:/Windows/Fonts/malgun.ttf")
  } else if (file.exists("/System/Library/Fonts/Supplemental/AppleGothic.ttf")) {
    font_add("nanum", "/System/Library/Fonts/Supplemental/AppleGothic.ttf")
  }
})

showtext_auto(TRUE)
options(ggtree_family = "nanum") 

## ===================== 1) phy 생성 & 라벨 매핑 =====================
phy <- as.phylo(fit)
map_label <- setNames(labs_mapped, labs_orig)

tmp <- unname(map_label[phy$tip.label])
phy$tip.label <- ifelse(is.na(tmp), sub("^X", "", phy$tip.label), tmp)

## 그룹 만들기
grp_raw <- cutree(fit, k = k)
names(grp_raw) <- sub("^X","", names(grp_raw))
names(grp_raw) <- ifelse(
  names(grp_raw) %in% names(map_P_to_name),
  unname(map_P_to_name[names(grp_raw)]),
  names(grp_raw)
)

grp_levels <- sort(unique(grp_raw))
tips_by_group <- lapply(grp_levels, function(g) names(grp_raw)[grp_raw == g])
names(tips_by_group) <- grp_levels

mrca_nodes <- sapply(tips_by_group, function(tips) if (length(tips) < 2) NA_integer_ else getMRCA(phy, tips))

## ===================== 2) 팔레트/스타일 =====================
fill_cols <- c("#7DB9B6", "#FDAF75", "#A3C4F3", "#C3DFA6")
if (length(grp_levels) != length(fill_cols)) {
  fill_cols <- scales::hue_pal()(length(grp_levels))
}
names(fill_cols) <- grp_levels

ring_outer_pad <- 0.6   
ring_thick     <- 1.0   

label_size   <- 2.6
label_offset <- 0.9

## ===================== 3) 트리(선)부터: groupOTU + geom_tree =====================
phy_grouped <- groupOTU(phy, tips_by_group)
p_circ <- ggtree(phy_grouped, layout = "circular", aes(color = group)) +
  geom_tree(size = 0.35) +
  scale_color_manual(values = fill_cols, na.value = "transparent") +
  # CHANGED: theme_void() completely removes the gray background grid lines
  theme_void(base_family = getOption("ggtree_family")) + 
  theme(legend.position = "none", plot.margin = margin(8,8,8,8))

## ===================== 4) 하이라이트 =====================
hl_formals <- names(formals(ggtree::geom_hilight))
support_offset <- "offset" %in% hl_formals

if (support_offset) {
  hl_offset <- 0.5
  hl_extend <- ring_thick
  for (g in names(mrca_nodes)) {
    nd <- mrca_nodes[[g]]; if (is.na(nd)) next
    p_circ <- p_circ +
      geom_hilight(node = nd, fill = fill_cols[g], colour = fill_cols[g],
                   alpha = 0.35, offset = hl_offset, extend = hl_extend)
  }
} else {
  d_tree <- p_circ$data
  R_tip  <- max(d_tree$x, na.rm = TRUE)
  R_out  <- R_tip + ring_outer_pad
  for (g in names(mrca_nodes)) {
    nd <- mrca_nodes[[g]]; if (is.na(nd)) next
    pos <- get_clade_position(p_circ, nd)  
    p_circ <- p_circ +
      annotate("rect", xmin = R_out - ring_thick, xmax = R_out,
               ymin = pos$ymin, ymax = pos$ymax,
               fill = fill_cols[g], color = fill_cols[g], alpha = 0.35)
  }
  p_circ <- p_circ + xlim(0, R_out + label_offset + 0.5)
}

## ===================== 5) 라벨: 기본은 검정, 한글 포함 시 빨강 =====================
# 1. Base labels (Default to black)
p_circ <- p_circ +
  geom_tiplab(aes(label = label),
              size = label_size, offset = label_offset,
              color = "black", family = getOption("ggtree_family"))

tip_df  <- p_circ$data %>% filter(isTip)

# 2. CHANGED: Find ALL labels that contain Korean characters [가-힣]
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label)) 

# 3. CHANGED: Overwrite Korean labels in RED
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# CHANGED: showtext_opts ensures vector conversion, removing cairo_pdf prevents conflict
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_2.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_2.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)
	   
=========================== 
==========================

## ===================== 0) 패키지 & 폰트 =====================
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2)
  library(ggtree); library(showtext); library(sysfonts)
})

# The most robust method for PDFs: Download Nanum Gothic directly from Google Fonts.
tryCatch({
  font_add_google("Nanum Gothic", "nanum")
}, error = function(e) {
  # Fallback if no internet
  if (file.exists("C:/Windows/Fonts/malgun.ttf")) {
    font_add("nanum", "C:/Windows/Fonts/malgun.ttf")
  } else if (file.exists("/System/Library/Fonts/Supplemental/AppleGothic.ttf")) {
    font_add("nanum", "/System/Library/Fonts/Supplemental/AppleGothic.ttf")
  }
})

showtext_auto(TRUE)
options(ggtree_family = "nanum") 

## ===================== 1) phy 생성 & 라벨 매핑 =====================
phy <- as.phylo(fit)
map_label <- setNames(labs_mapped, labs_orig)

tmp <- unname(map_label[phy$tip.label])
phy$tip.label <- ifelse(is.na(tmp), sub("^X", "", phy$tip.label), tmp)

## 그룹 만들기
grp_raw <- cutree(fit, k = k)
names(grp_raw) <- sub("^X","", names(grp_raw))
names(grp_raw) <- ifelse(
  names(grp_raw) %in% names(map_P_to_name),
  unname(map_P_to_name[names(grp_raw)]),
  names(grp_raw)
)

grp_levels <- sort(unique(grp_raw))
tips_by_group <- lapply(grp_levels, function(g) names(grp_raw)[grp_raw == g])
names(tips_by_group) <- grp_levels

mrca_nodes <- sapply(tips_by_group, function(tips) if (length(tips) < 2) NA_integer_ else getMRCA(phy, tips))

## ===================== 2) 팔레트/스타일 =====================
fill_cols <- c("#7DB9B6", "#FDAF75", "#A3C4F3", "#C3DFA6")
if (length(grp_levels) != length(fill_cols)) {
  fill_cols <- scales::hue_pal()(length(grp_levels))
}
names(fill_cols) <- grp_levels

ring_outer_pad <- 0.6   
ring_thick     <- 1.0   

label_size   <- 2.6
label_offset <- 0.9



## ===================== 3) 트리(선)부터: groupOTU + geom_tree =====================
phy_grouped <- groupOTU(phy, tips_by_group)

p_circ <- ggtree(phy_grouped, layout = "circular") +
  
  # 1. THE FIX: Draw a solid base tree first to bridge any coordinate gaps
  geom_tree(color = "gray80", linewidth = 0.6) + 
  
  # 2. Draw your colored grouped tree directly on top of it
  geom_tree(aes(color = group), linewidth = 0.6, lineend = "round", linejoin = "round") +
  
  scale_color_manual(values = fill_cols, na.value = "transparent") +
  theme_void(base_family = getOption("ggtree_family")) +
  theme(
    legend.position = "none", 
    plot.margin = margin(8,8,8,8),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

## ===================== 4) 하이라이트 (이전 코드와 동일하게 유지) =====================
# (Keep your existing Section 4 code here)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final4.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

========================= 

## ===================== 4) 하이라이트 (그룹 색상과 동일한 배경) =====================
# Loop through each group's Most Recent Common Ancestor (MRCA) node
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  # Only draw a highlight if the group has a valid ancestor node (2 or more tips)
  if (!is.na(nd)) {
    p_circ <- p_circ +
      geom_hilight(
        node = nd, 
        fill = fill_cols[g],   # Matches the highlight background to the branch color
        alpha = 0.2,           # Transparency (0.0 to 1.0) so lines/text aren't hidden
        extend = 0.5           # Pushes the highlight slightly past the tips for a clean border
      )
  }
}


## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final4.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

============================== 
## ===================== 4) 하이라이트 (균일한 바깥 링 형태) =====================
# 1. Calculate the maximum radius (the very edge of the tree tips)
d_tree <- p_circ$data
R_tip  <- max(d_tree$x, na.rm = TRUE)
R_out  <- R_tip + ring_outer_pad

# 2. Draw a uniform ring for each group
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    # Get the angular limits (ymin, ymax) of the specific clade
    pos <- get_clade_position(p_circ, nd)  
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = R_out - ring_thick, xmax = R_out, # Fixed width for all groups
               ymin = pos$ymin, ymax = pos$ymax,        # Angle spanned by the group
               fill = fill_cols[g], 
               color = "transparent",                   # No harsh borders
               alpha = 0.2)
  }
}

# 3. Expand the plot boundary so the labels don't get cut off by the new ring
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)


## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final5.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

============================ 
## ===================== 4) 하이라이트 (시작점과 끝점이 모두 동일한 웻지) =====================
d_tree <- p_circ$data

# 1. Find the starting point (inner level) and ending point (outer level)
# Get the x-coordinates of all group ancestor nodes
mrca_x <- sapply(mrca_nodes, function(nd) {
  if (is.na(nd)) return(NA)
  d_tree$x[d_tree$node == nd]
})

# Set a uniform inner level (starts at the deepest group's node)
# *If you want the inner ring closer to the center, you can manually set this (e.g., fixed_xmin <- 5)*
fixed_xmin <- min(mrca_x, na.rm = TRUE) 

# Set a uniform outer level (slightly past the tips)
R_tip <- max(d_tree$x, na.rm = TRUE)
R_out <- R_tip + 0.5 

# 2. Draw the uniform highlight wedges
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    # Get the angular limits (ymin, ymax) of the specific group
    pos <- get_clade_position(p_circ, nd)  
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = fixed_xmin, ymax = pos$ymax,  # All groups start at fixed_xmin
               xmax = R_out,      ymin = pos$ymin,  # All groups end at R_out
               fill = fill_cols[g], 
               color = "transparent", 
               alpha = 0.2)
  }
}

# 3. Expand the plot boundary so the labels fit nicely outside the new highlights
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final6.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	


============================ 
## ===================== 4) 하이라이트 (단일 균일 웻지) =====================
# (Make sure all previous geom_hilight code is completely deleted!)

d_tree <- p_circ$data

# 1. Calculate the single inner start line and outer end line
mrca_x <- sapply(mrca_nodes, function(nd) {
  if (is.na(nd)) return(NA)
  d_tree$x[d_tree$node == nd]
})

fixed_xmin <- min(mrca_x, na.rm = TRUE) # The uniform inner starting point
R_tip <- max(d_tree$x, na.rm = TRUE)
R_out <- R_tip + 0.5                    # The uniform outer ending point

# 2. Draw exactly ONE highlight wedge per group
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    pos <- get_clade_position(p_circ, nd)  
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = fixed_xmin, xmax = R_out,     # Every group uses these exact same X limits
               ymin = pos$ymin,   ymax = pos$ymax,  # The angle for the specific group
               fill = fill_cols[g], 
               color = "transparent",               # No border lines
               alpha = 0.2)                         # Keeps the single layer soft and readable
  }
}

# 3. Adjust the plot boundaries
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final7.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

=========================== 
## ===================== 4) 하이라이트 (시작점: 그룹 노드, 끝점: 균일한 바깥 링) =====================
d_tree <- p_circ$data

# 1. Calculate the uniform outer ending point (a perfect circle outside)
R_tip <- max(d_tree$x, na.rm = TRUE)
R_out <- R_tip + 0.5 

# 2. Draw the highlight wedges
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    # Get the angle (ymin, ymax) of the clade
    pos <- get_clade_position(p_circ, nd)  
    
    # Get the exact starting coordinate (x) of this specific group's root node
    node_x <- d_tree$x[d_tree$node == nd] 
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = node_x,     # STARTS at the group's specific root node
               xmax = R_out,      # ENDS at the uniform outer boundary
               ymin = pos$ymin,   ymax = pos$ymax,  
               fill = fill_cols[g], 
               color = "transparent",               
               alpha = 0.35)      # Richer color to match the "darker" look
  }
}

# 3. Adjust the plot boundaries so labels fit
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final8.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

=========================== 
## ===================== 4) 하이라이트 (완벽한 도넛 형태 웻지) =====================
d_tree <- p_circ$data

# 1. Set the perfect circular inner and outer boundaries
R_tip <- max(d_tree$x, na.rm = TRUE)
R_out <- R_tip + 0.5              # The uniform outer boundary

# THIS IS YOUR BLACK LINE: Adjust this number (e.g., 0.3, 0.4, 0.5) to change the inner circle size
R_in <- R_tip * 0.3               

# 2. Draw the continuous highlight ring
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    # Get the exact angle (ymin, ymax) for the group so they touch perfectly
    pos <- get_clade_position(p_circ, nd)  
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = R_in,       # All colors START at the exact same inner circle
               xmax = R_out,      # All colors END at the exact same outer circle
               ymin = pos$ymin,   ymax = pos$ymax,  
               fill = fill_cols[g], 
               color = "transparent",               
               alpha = 0.25)      # Soft pastel transparency
  }
}

# 3. Adjust plot boundaries to fit the labels outside the ring
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. English/Numeric ONLY (No Korean characters) -> Black
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = label_size, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. Korean ONLY -> Red and Bold
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = label_size + 0.2, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
# showtext 옵션을 활성화하여 폰트를 벡터 이미지로 완벽히 구워냅니다.
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_highlight_circ_v3_Final9.pdf", p_circ,
       width = 12, height = 12)
ggsave("Figure4_dendrogram_highlight_circ_v3_Final.png", p_circ,
       width = 3000, height = 3000, units = "px", dpi = 300)	

============================ 
## ===================== 0) 패키지 & 폰트 =====================
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(ggplot2)
  library(ggtree); library(showtext); library(sysfonts)
})

tryCatch({
  font_add_google("Nanum Gothic", "nanum")
}, error = function(e) {
  if (file.exists("C:/Windows/Fonts/malgun.ttf")) font_add("nanum", "C:/Windows/Fonts/malgun.ttf")
})

showtext_auto(TRUE)
options(ggtree_family = "nanum") 

## ===================== 1) phy 생성 & 라벨 매핑 =====================
phy <- as.phylo(fit)
map_label <- setNames(labs_mapped, labs_orig)

tmp <- unname(map_label[phy$tip.label])
phy$tip.label <- ifelse(is.na(tmp), sub("^X", "", phy$tip.label), tmp)

grp_raw <- cutree(fit, k = k)
names(grp_raw) <- sub("^X","", names(grp_raw))
names(grp_raw) <- ifelse(
  names(grp_raw) %in% names(map_P_to_name),
  unname(map_P_to_name[names(grp_raw)]),
  names(grp_raw)
)

grp_levels <- sort(unique(grp_raw))
tips_by_group <- lapply(grp_levels, function(g) names(grp_raw)[grp_raw == g])
names(tips_by_group) <- grp_levels

mrca_nodes <- sapply(tips_by_group, function(tips) if (length(tips) < 2) NA_integer_ else getMRCA(phy, tips))

## ===================== 2) 팔레트/스타일 =====================
fill_cols <- c("#7DB9B6", "#FDAF75", "#A3C4F3", "#C3DFA6")
if (length(grp_levels) != length(fill_cols)) {
  fill_cols <- scales::hue_pal()(length(grp_levels))
}
names(fill_cols) <- grp_levels

## ===================== 3) 트리(선) 그리기 =====================
phy_grouped <- groupOTU(phy, tips_by_group)

p_circ <- ggtree(phy_grouped, layout = "circular") +
  # Base gray tree
  geom_tree(color = "gray80", linewidth = 0.5) + 
  # Colored branches on top
  geom_tree(aes(color = group), linewidth = 0.5, lineend = "round", linejoin = "round") +
  
  scale_color_manual(values = fill_cols, na.value = "transparent") +
  theme_void(base_family = getOption("ggtree_family")) +
  theme(legend.position = "none", plot.margin = margin(8,8,8,8))

## ===================== 4) 하이라이트 (완벽한 도넛 & 검은선) =====================
d_tree <- p_circ$data
R_tip <- max(d_tree$x[d_tree$isTip == TRUE], na.rm = TRUE)

# -- 치수 설정 (Adjust these if you want the ring thicker/thinner) --
R_in  <- R_tip * 0.45   # Inner circle (The black line position)
R_out <- R_tip + 0.8    # Outer circle (Where the highlight ends)
label_offset <- 1.0     # Pushes the text safely outside the color ring

# 1. 완벽하게 이어지는 도넛 배경 그리기
for (g in names(mrca_nodes)) {
  nd <- mrca_nodes[[g]]
  
  if (!is.na(nd)) {
    pos <- get_clade_position(p_circ, nd)  
    
    p_circ <- p_circ +
      annotate("rect", 
               xmin = R_in, xmax = R_out, 
               # -0.5 and +0.5 mathematically closes the white gaps between groups!
               ymin = pos$ymin - 0.5, ymax = pos$ymax + 0.5,  
               fill = fill_cols[g], 
               color = "transparent",               
               alpha = 0.25)
  }
}

# 2. 안쪽 검은색 원 그리기 (The Black Circle from your image)
# In polar coordinates, a vertical line (xintercept) becomes a perfect circle!
p_circ <- p_circ + geom_vline(xintercept = R_in, color = "black", linewidth = 0.6)

# 확장된 여백 적용
p_circ <- p_circ + xlim(0, R_out + label_offset + 1.0)

## ===================== 5) 라벨: 완벽 분리 (영어=검정, 한글=빨강) =====================
tip_df <- p_circ$data %>% filter(isTip)

# 1. 영어/숫자 (검정)
tip_english <- tip_df %>% filter(!grepl("[가-힣]", label))
if (nrow(tip_english) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_english, aes(label = label),
                size = 2.6, offset = label_offset,
                color = "black", family = getOption("ggtree_family"))
}

# 2. 한글 (빨강 & 굵게)
tip_korean <- tip_df %>% filter(grepl("[가-힣]", label))
if (nrow(tip_korean) > 0) {
  p_circ <- p_circ +
    geom_tiplab(data = tip_korean, aes(label = label),
                size = 2.8, offset = label_offset,
                color = "red", fontface = "bold",
                family = getOption("ggtree_family"))
}

## ===================== 6) 저장 =====================
showtext_opts(dpi = 300) 

ggsave("Figure4_dendrogram_PerfectDonut.pdf", p_circ, width = 12, height = 12)
ggsave("Figure4_dendrogram_PerfectDonut.png", p_circ, width = 3000, height = 3000, units = "px", dpi = 300)


==========Fst========= 
# Clear the file if it already exists
> pop_assignments.txt

# Process each group and append to the final file
awk '{print $1, $2, "GroupA"}' groupA.fam >> pop_assignments.txt
awk '{print $1, $2, "GroupB"}' groupB.fam >> pop_assignments.txt
awk '{print $1, $2, "GroupC"}' groupC.fam >> pop_assignments.txt
awk '{print $1, $2, "GroupD"}' groupD.fam >> pop_assignments.txt

plink --bfile input_data.plink --within pop_assignments.txt --fst --out final_fst_results

# Sort by the Fst column (column 5) in descending order and show the top 10
sort -k5 -gr final_fst_results.fst | head -n 11

#Calculate Pairwise Fst
#conda install -c bioconda plink2
#plink2 --version

# Create a new file with a header (FID, IID, and POP)
echo -e "FID\tIID\tPOP" > pop_assignments_header.txt
cat pop_assignments.txt >> pop_assignments_header.txt

plink2 --bfile input_data.plink \
       --pheno pop_assignments_header.txt \
       --fst POP \
       --out pairwise_results
	   
#Create the Matrix in R
# Load libraries
library(ggplot2)
library(reshape2)

# Create the data frame from your results
fst_data <- data.frame(
  P1 = c("GroupA", "GroupA", "GroupA", "GroupB", "GroupB", "GroupC"),
  P2 = c("GroupB", "GroupC", "GroupD", "GroupC", "GroupD", "GroupD"),
  Fst = c(0.391216, 0.362373, 0.681005, 0.363502, 0.674467, 0.697982)
)

# Get all unique group names
groups <- unique(c(fst_data$P1, fst_data$P2))
n <- length(groups)

# Initialize an empty matrix
fst_matrix <- matrix(0, nrow=n, ncol=n, dimnames=list(groups, groups))

# Fill the matrix (both upper and lower triangles)
for(i in 1:nrow(fst_data)) {
  fst_matrix[fst_data$P1[i], fst_data$P2[i]] <- fst_data$Fst[i]
  fst_matrix[fst_data$P2[i], fst_data$P1[i]] <- fst_data$Fst[i]
}

# Print the matrix to console
print(fst_matrix)

# Save the matrix to a CSV file
write.csv(fst_matrix, "pairwise_fst_matrix.csv")

#Visualize with a Heatmap
# Melt the matrix for ggplot
melted_fst <- melt(fst_matrix)

# Create heatmap (plot_fst.R)
# Load libraries
library(ggplot2)
library(reshape2)

# 1. Prepare the data
fst_data <- data.frame(
  P1 = c("GroupA", "GroupA", "GroupA", "GroupB", "GroupB", "GroupC"),
  P2 = c("GroupB", "GroupC", "GroupD", "GroupC", "GroupD", "GroupD"),
  Fst = c(0.391216, 0.362373, 0.681005, 0.363502, 0.674467, 0.697982)
)

groups <- unique(c(fst_data$P1, fst_data$P2))
n <- length(groups)
fst_matrix <- matrix(0, nrow=n, ncol=n, dimnames=list(groups, groups))

for(i in 1:nrow(fst_data)) {
  fst_matrix[fst_data$P1[i], fst_data$P2[i]] <- fst_data$Fst[i]
  fst_matrix[fst_data$P2[i], fst_data$P1[i]] <- fst_data$Fst[i]
}

melted_fst <- melt(fst_matrix)

# 2. Create the plot object
p <- ggplot(data = melted_fst, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#fff5f5", high = "#ff0000") +
  geom_text(aes(label = sprintf("%.3f", value)), size = 5) +
  theme_minimal() +
  labs(title = "Pairwise Fst Genetic Distance", 
       subtitle = "Group D shows highest divergence",
       x = "", y = "", fill = "Fst Value") +
  theme(axis.text = element_text(size = 12, face = "bold"))

# 3. Save the file
# Save as PNG for presentations
ggsave("fst_heatmap.png", plot = p, width = 8, height = 6, dpi = 300)

# Save as PDF for publications (vector format)
ggsave("fst_heatmap.pdf", plot = p, width = 8, height = 6)

cat("Success: Heatmap saved as fst_heatmap.png and fst_heatmap.pdf\n")





plink2 --bfile input_data.plink \
       --pheno pop_assignments.txt \
       --fst \
       --out pairwise_results


#plot_manhattan.R	   
# Load the library
if (!require("qqman")) {
    install.packages("qqman", repos="http://cran.us.r-project.org")
    library(qqman)
}

# Load the data from your previous PLINK run
data <- read.table("final_fst_results.fst", header=TRUE)

# Remove markers with missing Fst values
data <- na.omit(data)

# Open a PNG file device to save the plot
# You can change 'res' for higher quality or 'width'/'height' for size
png("fst_manhattan_plot.png", width=1200, height=600, res=120)

# Generate the Manhattan plot
# We use logp=FALSE because Fst values are not p-values
manhattan(data, 
          chr="CHR", 
          bp="POS", 
          snp="SNP", 
          p="FST", 
          logp=FALSE, 
          ylab="Fst Value", 
          main="Global Fst Across the Genome",
          col=c("blue4", "orange3"), # Alternating chromosome colors
          suggestiveline=FALSE,      # Removes the horizontal p-value lines
          genomewideline=FALSE)

# Close the device to finish saving the file
dev.off()

cat("Manhattan plot saved as: fst_manhattan_plot.png\n")


# Get the top 10 SNPs with the highest Fst values
sort -k5 -gr final_fst_results.fst | head -n 11


(sra_env) [krjo@dongph_genome cc]$ sort -k5 -gr final_fst_results.fst | head -n 11
7	Chr7_9941999	9941999	231	1
7	Chr7_9803370	9803370	231	1
7	Chr7_9722716	9722716	231	1
7	Chr7_12682956	12682956	231	1
7	Chr7_11042802	11042802	231	1
7	Chr7_10826629	10826629	231	1
7	Chr7_10658542	10658542	231	1
7	Chr7_10012470	10012470	231	1
6	Chr6_22933739	22933739	231	1
6	Chr6_17443569	17443569	231	1
6	Chr6_16349403	16349403	231	1
(sra_env) [krjo@dongph_genome cc]$ 
#######################################
#This output shows your **Top Outlier SNPs**—the genetic markers that are most different between your four groups.
#
#In your data, these SNPs have an ** value of 1.0**. This is a very strong result.
#
### 1. What an  of 1.0 means
#
#An  of 1.0 means **complete fixation**. At these specific locations in the genome:

* One group (or set of groups) has one version of the DNA (e.g., all individuals have "A").
* Another group has a completely different version (e.g., all individuals have "G").
* There are **zero** shared genotypes between them at these positions.

### 2. Genetic "Hotspots" on Chromosome 7 and 6

Looking at your results, you have found major "selection peaks" or "islands of differentiation":

* **Chromosome 7 Cluster:** Notice that several top SNPs (Chr7_9941999, Chr7_9803370, etc.) are physically very close to each other (around the 9MB to 12MB region). This suggests there is a **large block of DNA on Chromosome 7** that is completely different between your populations.
* **Chromosome 6 Cluster:** Similarly, you have markers around 16MB to 22MB on Chromosome 6 that show complete differentiation.

### 3. Biological Interpretation

When you see a cluster of SNPs with  like this, it usually indicates one of two things:

1. **Strong Natural Selection:** A gene in this region is so important for the survival of one group (e.g., adapting to a specific environment, soil, or climate) that 100% of the population carries the beneficial version.
2. **Speciation/Long Isolation:** The groups have been separated for so long without interbreeding that these parts of the genome have drifted entirely apart. Given your weighted  of 0.55, these groups are behaving more like distinct subspecies or varieties than a single interbreeding population.

### Next Step: Identifying the Genes

Since these SNPs are "fixed" (1.0), they are likely linked to the traits that define your groups. To find out what these traits are:

1. Go to a genome browser for your species (e.g., NCBI or a specific crop database).
2. Search for the region **Chr7: 9,000,000 - 13,000,000**.
3. Look for genes in that window.

**Would you like me to help you write a script to extract the exact genotypes at these "fixed" positions so you can see which group has which allele?**
#
#
######################################
# Basic Fst calculation
plink --bfile input_data.plink \
      --within pop_assignments.txt \
      --fst \
      --out fst_results



# Calculate Fst for all population pairs
plink2 --bfile input_data.plink \
       --pheno pop_assignments.txt \
       --fst \
       --out fst_results_p2



# 1. Convert PLINK to VCF
plink --bfile input_data.plink --recode vcf --out temp_data

# 2. Create population list files (one ID per line)
# pop1.txt contains IDs for group 1
# pop2.txt contains IDs for group 2

# 3. Calculate Fst in 50kb windows
vcftools --vcf temp_data.vcf \
         --weir-fst-pop pop1.txt \
         --weir-fst-pop pop2.txt \
         --fst-window-size 50000 \
         --fst-window-step 10000 \
         --out windowed_fst
		 
=======plot_manhattan_highlight.R========== 
# Load library
if (!require("qqman")) {
    install.packages("qqman", repos="http://cran.us.r-project.org")
    library(qqman)
}

# Load the data
data <- read.table("final_fst_results.fst", header=TRUE)
data <- na.omit(data)

# 1. Define the list of SNPs to highlight based on your results
snps_to_highlight <- c(
  "Chr7_9941999", "Chr7_9803370", "Chr7_9722716", "Chr7_12682956", 
  "Chr7_11042802", "Chr7_10826629", "Chr7_10658542", "Chr7_10012470",
  "Chr6_22933739", "Chr6_17443569", "Chr6_16349403"
)

# 2. Open PNG device
png("fst_manhattan_highlighted.png", width=1200, height=600, res=120)

# 3. Create the Manhattan plot with highlighting
manhattan(data, 
          chr="CHR", 
          bp="POS", 
          snp="SNP", 
          p="FST", 
          logp=FALSE, 
          highlight = snps_to_highlight, # This circles the specific SNPs
          col=c("grey60", "grey30"),     # Muted background colors
          suggestiveline=0.8,            # Add a dashed line at Fst 0.8
          genomewideline=FALSE,
          ylab="Fst Value", 
          main="Fst Across Genome (Fixed Outliers Highlighted)")

# 4. Close device
dev.off()

cat("Highlighted Manhattan plot saved as: fst_manhattan_highlighted.png\n")

=========================== 
#extract the exact genotypes at these highly differentiated (fixed) positions
#Create a list of the target SNPs
cat << EOF > fixed_snps.txt
Chr2_245336
Chr2_2425665
Chr2_2710369
Chr2_3050360
Chr2_3846970
Chr2_3924266
Chr2_4525854
Chr2_4608925
Chr2_5654644
Chr2_6045153
Chr2_6395179
Chr2_18558011
Chr3_2462884
Chr3_2566140
Chr3_2942055
Chr3_3137604
Chr3_4910770
Chr3_5254526
Chr3_5304137
Chr3_21209622
Chr4_8835511
Chr4_9254036
Chr4_12291355
Chr4_12944658
Chr4_20594356
Chr4_21646544
Chr4_21826474
Chr6_16163127
Chr6_16349403
Chr6_17443569
Chr6_22933739
Chr7_9722716
Chr7_9803370
Chr7_9941999
Chr7_10012470
Chr7_10658542
Chr7_10826629
Chr7_11042802
Chr7_12682956
EOF

#Extract Genotypes using PLINK
# Extract these SNPs into a readable format
plink --bfile input_data.plink \
      --extract fixed_snps.txt \
      --recode \
      --out fixed_genotypes
	  
#Summarize by Group (summarize_fixed.py)	  
import pandas as pd

# Load assignments
pops = pd.read_csv('pop_assignments.txt', sep='\s+', header=None, names=['FID', 'IID', 'Group'])

# Load the map file to get SNP order
names = ['CHR', 'SNP', 'DIST', 'POS']
map_df = pd.read_csv('fixed_genotypes.map', sep='\s+', header=None, names=names)

# Load the ped file (Genotypes start at column 7)
ped = pd.read_csv('fixed_genotypes.ped', sep='\s+', header=None)

# Extract just IDs and Genotypes
genotypes = ped.iloc[:, [1] + list(range(6, ped.shape[1]))]
genotypes.columns = ['IID'] + [f"{s}_{i}" for s in map_df['SNP'] for i in [1, 2]]

# Merge with population info
final = pd.merge(pops[['IID', 'Group']], genotypes, on='IID')

# Calculate the most frequent genotype per group for the first few SNPs
for snp in map_df['SNP'][:5]:
    print(f"\n--- SNP: {snp} ---")
    summary = final.groupby('Group')[[f"{snp}_1", f"{snp}_2"]].first()
    print(summary)
	
	
#create a comprehensive summary of your Fst analysis
# 1. Create the header
echo "===============================================" > fst_summary_report.txt
echo "       GENETIC DIFFERENTIATION SUMMARY        " >> fst_summary_report.txt
echo "===============================================" >> fst_summary_report.txt
echo "" >> fst_summary_report.txt

# 2. Extract Global Fst Statistics from the log file
echo "[1] GLOBAL FST STATISTICS" >> fst_summary_report.txt
grep -E "Mean Fst estimate|Weighted Fst estimate" final_fst_results.log >> fst_summary_report.txt
echo "" >> fst_summary_report.txt

# 3. Add Pairwise Fst Table
echo "[2] PAIRWISE FST (HUDSON METHOD)" >> fst_summary_report.txt
if [ -f "pairwise_results.fst.summary" ]; then
    cat pairwise_results.fst.summary >> fst_summary_report.txt
else
    echo "Pairwise results file not found." >> fst_summary_report.txt
fi
echo "" >> fst_summary_report.txt

# 4. Extract Top 10 Outlier SNPs (Fixed Differences)
echo "[3] TOP 10 OUTLIER SNPs (FIXED DIFFERENCES)" >> fst_summary_report.txt
echo -e "CHR\tSNP\tPOS\tNMISS\tFST" >> fst_summary_report.txt
sort -k5 -gr final_fst_results.fst | head -n 10 >> fst_summary_report.txt
echo "" >> fst_summary_report.txt

# 5. Summary of Significant Clusters
echo "[4] GENOMIC HOTSPOTS IDENTIFIED" >> fst_summary_report.txt
echo " - Chromosome 7: Large cluster of fixed SNPs (~9MB - 12MB)" >> fst_summary_report.txt
echo " - Chromosome 6: Significant differentiation observed (~16MB - 22MB)" >> fst_summary_report.txt
echo "" >> fst_summary_report.txt

echo "Report generated: fst_summary_report.txt"

================= 
# 1. Start with Global Statistics
echo "SECTION,METRIC,VALUE" > fst_comprehensive_summary.csv
echo "Global,Mean Fst,0.501984" >> fst_comprehensive_summary.csv
echo "Global,Weighted Fst,0.550057" >> fst_comprehensive_summary.csv
echo "Global,Valid Markers,998" >> fst_comprehensive_summary.csv
echo "Global,Excluded Markers,104" >> fst_comprehensive_summary.csv
echo ",," >> fst_comprehensive_summary.csv

# 2. Add Pairwise Fst Comparison
echo "SECTION,POP1,POP2,HUDSON_FST" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupA,GroupB,0.391216" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupA,GroupC,0.362373" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupA,GroupD,0.681005" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupB,GroupC,0.363502" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupB,GroupD,0.674467" >> fst_comprehensive_summary.csv
echo "Pairwise,GroupC,GroupD,0.697982" >> fst_comprehensive_summary.csv
echo ",,," >> fst_comprehensive_summary.csv

# 3. Add Top 10 Outlier SNPs (The "Fixed" SNPs)
echo "SECTION,CHR,SNP,POS,FST_VALUE" >> fst_comprehensive_summary.csv
# This pulls directly from your results file and formats it for CSV
sort -k5 -gr final_fst_results.fst | head -n 10 | awk '{print "Outlier,"$1","$2","$3","$5}' >> fst_comprehensive_summary.csv

echo "Excel-friendly summary created: fst_comprehensive_summary.csv"	

###extract all SNPs where the $F_{st}$ value is exactly 1.0
# Filter for Fst = 1.0 (Column 5) and save as CSV
awk 'BEGIN {OFS=","} $5 == 1.0 {print $1, $2, $3, $5}' final_fst_results.fst > fixed_snps_list.csv

# Just extract the SNP names (Column 2)
awk '$5 == 1.0 {print $2}' final_fst_results.fst > fixed_snp_names.txt

# Use the list to create a new small PLINK dataset of only fixed SNPs
plink --bfile input_data.plink --extract fixed_snp_names.txt --make-bed --out fixed_only_data

echo "Total number of fixed SNPs (Fst=1.0):"
awk '$5 == 1.0' final_fst_results.fst | wc -l

awk '$5 == 1.0 {print $2}' final_fst_results.fst > truly_fixed_snps.txt

#================================= 
#Step 1: Calculate LD in PLINK
# Calculate LD for groupA
plink --bfile groupA \
      --ld-window-r2 0 \
      --ld-window 999999 \
      --ld-window-kb 1000 \
      --r2 gz \
      --out groupA_LD

#Step 2: Process and Plot in R
library(ggplot2)
library(dplyr)
library(data.table)

# Function to read and bin LD data
process_ld <- function(filename, group_name) {
  df <- fread(filename)
  # Calculate distance
  df$dist <- abs(df$BP_B - df$BP_A)
  # Create bins of 10kb
  df$dist_bin <- round(df$dist / 10000) * 10000
  
  # Summarize mean R2 per bin
  df %>%
    group_by(dist_bin) %>%
    summarise(mean_r2 = mean(R2)) %>%
    mutate(Group = group_name)
}

# Load all groups
dataA <- process_ld("groupA_LD.ld.gz", "Group A")
dataB <- process_ld("groupB_LD.ld.gz", "Group B")
dataC <- process_ld("groupC_LD.ld.gz", "Group C")
dataD <- process_ld("groupD_LD.ld.gz", "Group D")

# Combine and Plot
all_data <- rbind(dataA, dataB, dataC, dataD)

ggplot(all_data, aes(x = dist_bin / 1000, y = mean_r2, color = Group)) +
  geom_line(size = 1) +
  labs(x = "Distance (kb)", y = expression(r^2), title = "LD Decay Comparison") +
  theme_minimal()


#Extract Genotypes for Each Group
# Create a list of your group names
groups=("groupA" "groupB" "groupC" "groupD")

for grp in "${groups[@]}"; do
    # Extract the genotypes for this specific group
    # We use --keep to select only the IDs present in your .fam file
    plink --bfile input_data.plink \
          --keep ${grp}.fam \
          --make-bed \
          --out ${grp}_genotypes

    # Immediately calculate the LD decay for this group
    # This saves disk space by not keeping massive intermediate files
    plink --bfile ${grp}_genotypes \
          --ld-window-r2 0 \
          --ld-window 999999 \
          --ld-window-kb 1000 \
          --r2 gz \
          --out ${grp}_LD_results
done


#plink --bfile ${grp}_genotypes --maf 0.05 --ld-window-r2 0 ... (rest of command)

#The R Visualization Script(Visualization.R)

#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)

# 1. Define the files and group names
files <- c("groupA_LD_results.ld.gz", "groupB_LD_results.ld.gz", 
           "groupC_LD_results.ld.gz", "groupD_LD_results.ld.gz")
group_names <- c("Group A", "Group B", "Group C", "Group D")

# 2. Function to process each file
process_ld_decay <- function(file, name) {
  cat("Processing:", name, "...\n")
  df <- fread(file)
  
  # Calculate distance between SNPs
  df[, dist := abs(BP_B - BP_A)]
  
  # Bin distances by 10kb (10,000 bp)
  df[, dist_bin := round(dist / 10000) * 10000]
  
  # Calculate mean R2 per bin
  decay <- df[, .(mean_r2 = mean(R2)), by = dist_bin]
  decay$Group <- name
  return(decay)
}

# 3. Combine data from all four groups
all_groups_decay <- rbindlist(lapply(1:4, function(i) process_ld_decay(files[i], group_names[i])))

# 4. Create the Plot
ggplot(all_groups_decay, aes(x = dist_bin / 1000, y = mean_r2, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(alpha = 0.5) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "LD Decay Comparison Across Groups",
    x = "Physical Distance (kb)",
    y = expression(Average ~ r^2),
    color = "Population"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# 5. Save the plot to a file
ggsave("LD_Decay_Comparison.png", width = 10, height = 7, dpi = 300)

cat("Plot saved as LD_Decay_Comparison.png\n")  
  
###########Visualization_v2.R################
ggplot(all_groups_decay, aes(x = dist_bin / 1000, y = mean_r2, color = Group)) +
  # This adds the smooth trend line
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) + 
  # This keeps the raw data as faint points in the background
  geom_point(alpha = 0.2) + 
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Smoothed LD Decay Comparison",
    x = "Physical Distance (kb)",
    y = expression(Average ~ r^2)
  ) +
  theme_minimal()

# 5. Save the plot to a file
ggsave("LD_Decay_Comparison2.png", width = 10, height = 7, dpi = 300)

cat("Plot saved as LD_Decay_Comparison2.png\n")  
  
  
#Optional: R Script for Half-Decay Values
library(dplyr)

# 1. Identify the 'Half-Max' for each group based on the decay curve
half_decay_table <- all_groups_decay %>%
  group_by(Group) %>%
  arrange(dist_bin) %>%
  mutate(
    max_r2 = max(mean_r2),
    min_r2 = min(mean_r2),
    half_way_point = (max_r2 + min_r2) / 2
  ) %>%
  # Find the first bin where mean_r2 drops below the half_way_point
  filter(mean_r2 <= half_way_point) %>%
  slice(1) %>% 
  select(Group, half_decay_kb = dist_bin) %>%
  mutate(half_decay_kb = half_decay_kb / 1000)

print("--- Half-Decay Distances (kb) ---")
print(half_decay_table)

# 2. Export the table to a CSV for your report
write.csv(half_decay_table, "LD_Half_Decay_Results.csv", row.names = FALSE)

###Final Check: Sample Sizes
wc -l groupA.fam groupB.fam groupC.fam groupD.fam


#Balance Groups
# Set the target size based on Group D
TARGET=31

# Randomly sample from A, B, and C
for grp in "groupA" "groupB" "groupC"; do
    shuf -n $TARGET ${grp}.fam > ${grp}_subsampled.fam
done

# Copy Group D as it is already at 31
cp groupD.fam groupD_subsampled.fam

echo "Subsampling complete. All groups now have 31 inds."



#Re-run the PLINK Analysis
groups=("groupA" "groupB" "groupC" "groupD")

for grp in "${groups[@]}"; do
    # Extract genotypes for the subsampled 31 inds
    plink --bfile input_data.plink \
          --keep ${grp}_subsampled.fam \
          --make-bed \
          --out ${grp}_sub_genotypes

    # Calculate LD
    plink --bfile ${grp}_sub_genotypes \
          --ld-window-r2 0 \
          --ld-window-kb 1000 \
          --r2 gz \
          --out ${grp}_sub_LD
done

#LD_decay_Visualization_v3.R
#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)

# 1. Define the files and group names for the SUBSAMPLED data
# Note: Ensure these filenames match the output of your PLINK loop exactly.
files <- c("groupA_sub_LD.ld.gz", "groupB_sub_LD.ld.gz", 
           "groupC_sub_LD.ld.gz", "groupD_sub_LD.ld.gz")
group_names <- c("Group A", "Group B", "Group C", "Group D")

# 2. Function to process each file
process_ld_decay <- function(file, name) {
  if(!file.exists(file)){
    stop(paste("File not found:", file))
  }
  cat("Processing Subsampled:", name, "...\n")
  df <- fread(file)
  
  # Calculate distance between SNPs
  df[, dist := abs(BP_B - BP_A)]
  
  # Bin distances by 10kb (10,000 bp)
  df[, dist_bin := round(dist / 10000) * 10000]
  
  # Calculate mean R2 per bin
  decay <- df[, .(mean_r2 = mean(R2)), by = dist_bin]
  decay$Group <- name
  return(decay)
}

# 3. Combine data from all four groups
all_groups_decay <- rbindlist(lapply(1:4, function(i) process_ld_decay(files[i], group_names[i])))

# 4. Create the Balanced Plot
# Using linewidth as per the previous warning message
ggplot(all_groups_decay, aes(x = dist_bin / 1000, y = mean_r2, color = Group)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) + 
  geom_point(alpha = 0.2) + 
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Balanced LD Decay Comparison (n=31 per group)",
    subtitle = "Subsampled to control for sample size bias",
    x = "Physical Distance (kb)",
    y = expression(Average ~ r^2)
  ) +
  theme_minimal()

# 5. Save the updated plot
ggsave("LD_Decay_Comparison_Subsampled.png", width = 10, height = 7, dpi = 300)

# 6. Calculate the new balanced Half-Decay table
half_decay_table <- all_groups_decay %>%
  group_by(Group) %>%
  arrange(dist_bin) %>%
  mutate(
    max_r2 = max(mean_r2),
    min_r2 = min(mean_r2),
    half_way_point = (max_r2 + min_r2) / 2
  ) %>%
  filter(mean_r2 <= half_way_point) %>%
  slice(1) %>% 
  select(Group, half_decay_kb = dist_bin) %>%
  mutate(half_decay_kb = half_decay_kb / 1000)

write.csv(half_decay_table, "LD_Half_Decay_Subsampled.csv", row.names = FALSE)

cat("Balanced Analysis Complete.\nPlot: LD_Decay_Comparison_Subsampled.png\nTable: LD_Half_Decay_Subsampled.csv\n")












