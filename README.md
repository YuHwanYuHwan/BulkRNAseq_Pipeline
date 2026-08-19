# BulkRNAseq_Pipeline

raw FASTQ → count matrix → TMM/CPM. 전처리 전용. DEG·GO 는 별도 저장소.

전처리는 **이 저장소 한 곳에서만** 한다. 프로젝트별로 복제하지 않는다.
연구 프로젝트는 `Output/<project>/<group>/` 만 가져다 쓴다.

---

## 1. Setup

```bash
git clone <repo> && cd BulkRNAseq_Pipeline
bash setup.sh --create-env      # conda 환경 생성 + 전체 점검
```

`setup.sh` 는 **무엇이 빠졌는지 먼저 알려준다.** 6시간짜리 STAR 가 3단계에서 죽고 나서
참조 게놈이 없었다는 걸 알게 되는 상황을 막는 게 목적이다. 몇 번을 돌려도 안전하다 —
`--create-env` 없이는 아무것도 설치하지 않고 상태만 보고한다.

점검 항목: config → conda 환경 → 도구 8종 + edgeR → 참조 게놈(종별 FASTA·GTF) → 디스크 여유
→ 로직 self-check.

### 소프트웨어

`setup.sh --create-env` 가 아래를 만든다. 수동으로 하려면:

```bash
conda create -n bulkrnaseq -c bioconda -c conda-forge \
    sra-tools fastqc cutadapt star htseq multiqc \
    bioconductor-edger r-base
```

| 도구 | 용도 |
|---|---|
| sra-tools | 공공 데이터 다운로드 (`prefetch`, `fasterq-dump`) |
| FastQC | raw read QC |
| cutadapt | adapter trimming |
| STAR | alignment + index build |
| HTSeq | gene-level count |
| MultiQC | QC 집계 |
| R + edgeR | TMM/CPM 정규화 |

### config

`setup.sh` 가 `config.sh.example` 를 `config.sh` 로 복사한다. PATH 에 도구가 다 있으면
전부 비워둬도 된다. 이미 설치해둔 FastQC 를 쓰려면 `FASTQC_BIN` 을, 랩 공용 참조 게놈이
있으면 `REF_ROOT` 를 지정한다.

### 참조 게놈 (직접 다운로드)

`reference_Genomes/<종명>/` 아래에 FASTA 와 GTF 를 둔다. 종명은 **학명**으로 — Ensembl
파일명과 맞춰야 스크립트가 파일을 찾는다.

```
reference_Genomes/
  Homo_sapiens/
    Homo_sapiens.GRCh38.dna.primary_assembly.fa
    Homo_sapiens.GRCh38.113.gtf
  Mus_musculus/
    Mus_musculus.GRCm39.dna.primary_assembly.fa
    Mus_musculus.GRCm39.115.gtf
```

[Ensembl FTP](https://ftp.ensembl.org/pub/) 에서 받는다.
```bash
# 예: human release 113
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip *.gz
```

**STAR 인덱스는 받지 않는다.** 필요한 overhang 이 없으면 파이프라인이 자동으로 빌드하고
`reference_Genomes/index_registry.tsv` 에 기록한 뒤 영구 보존한다. 다음부터는 재사용된다.
인덱스 하나에 디스크 약 30GB, 빌드에 RAM 32GB 이상이 필요하다.

## 2. 사용법

### 데이터 배치

```
rawData/<project>/<group>/
    group.conf
    SRR27743399_1.fastq.gz      ← 파일명 stem = sample_id
    SRR27743399_2.fastq.gz
    GSM7712347/                 ← 폴더 = sample_id, 안의 run 들은 merge 대상
        SRR27743401_1.fastq.gz
        ...
```

**project** 는 연구 단위, **group** 은 count matrix 를 뽑는 단위다.
Musculoskeletal_Map 프로젝트 아래 sarcopenia · cachexia 처럼 나뉜다.

공공 데이터면 다운로드 스크립트가 배치까지 해준다. 한 샘플에 run 이 여러 개면 자동으로 폴더로 묶인다.

```bash
bash Scripts/PublicData_download.sh rawData/MSK_Map/sarcopenia SRR27743399 SRR27743400
bash Scripts/PublicData_download.sh rawData/MSK_Map/sarcopenia srr_list.txt
```

### group.conf

사람이 쓰는 파일은 이것 하나다.

```bash
species       = Homo_sapiens     # 필수. reference_Genomes/ 의 폴더명과 일치
adapter_kit   = Illumina_TruSeq  # 생략 시 Illumina_universal
strandedness  =                  # 처음엔 비워둔다. stage 1 후 probe 결과 보고 기입
```

### 실행

```bash
bash Scripts/run_stage1.sh rawData/MSK_Map/sarcopenia     # 또는 sbatch
```
FastQC → trimming → STAR → strandedness probe. probe 출력을 보고 `group.conf` 의
`strandedness` 에 `no` / `yes` / `reverse` 중 하나를 적는다.

```bash
bash Scripts/run_stage2.sh rawData/MSK_Map/sarcopenia
```
HTSeq → count matrix → CPM → MultiQC + 리포트.

개별 단계만 돌려도 된다. 모든 스크립트가 `<group_dir>` 하나만 받는다.

### 결과

```
Output/<project>/<group>/
    <group>_count_matrix.tsv      raw counts (Gene × Samples)
    <group>_CPM.tsv               TMM-normalized CPM
    <group>_multiqc_report.html
    <group>_pipeline_report.txt   버전·파라미터·QC + Methods 초안 문단
```

`Processed/` 는 중간물(trimmed FASTQ, BAM)이라 **언제든 삭제해도 된다.**
raw 와 이 스크립트만 있으면 재생성된다.

---

## 3. 설계 규칙

- **샘플 식별은 디렉토리가 진실원**이다. 시트를 따로 쓰지 않는다.
  폴더 = 샘플(merge), 평면 파일 = 개별 샘플, `_1`/`_2` = paired.
- **strandedness 는 추론하지 않는다.** kit 이름이 예측하지 못한다 —
  `SureSelect Strand Specific` 인데 unstranded 인 데이터가 실재한다.
  probe 로 확인하고 사람이 적는다. 비어 있으면 stage 2 가 시작을 거부한다.
- **sjdbOverhang 은 데이터에서 산출**한다. trimmed read 의 최대 길이 − 1.
  group 전체가 같은 인덱스를 쓴다 — 샘플마다 다르면 count 를 비교할 수 없다.
- **`.done` 파일로 재개**한다. 중간에 죽어도 끝난 샘플은 건너뛴다.

## 4. 새 kit 추가

`Scripts/AdapterSequenceList.csv` 에 한 줄 붙이면 이후 모든 데이터에 재사용된다.

```csv
kit,adapter_R1,adapter_R2
MyNewKit,AGATCGGAAGAGCACACGTCT,AGATCGGAAGAGCGTCGTGTA
```

## 5. 자체 점검

```bash
bash setup.sh                      # 환경 전체 (self-check 포함)
for f in Scripts/tests/*.sh; do bash "$f"; done   # 로직만
```

샘플 스캔·merge 그룹핑·overhang 산출·matrix 조립·probe 판정 로직을 검증한다.
**생물정보 도구가 하나도 없어도 돈다** — clone 직후 자기 환경에서 로직이 정상인지
먼저 확인할 수 있다.
