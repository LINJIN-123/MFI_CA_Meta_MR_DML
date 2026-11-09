# 代码复现、永久归档与可引用 DOI 操作指南（科研新手版）

> 适用论文：含 meta 分析（R）、双向两样本孟德尔随机化 2SMR（R）、汇总水平双机器学习 DML（Python）的流行病学论文。
> 本指南配套的本地仓库已经按规范建好，可直接使用。所有命令均给出图形化（优先）与命令行（备选）两种方式。
> 核心原则：**GitHub 是"开发仓库"，Zenodo 才是"带 DOI 的正式永久归档"，两者分工不同，不能互相替代。**

---

## 模块 1　本地仓库目录结构与环境锁定

### 1.1 仓库目录结构（已建好，根目录 `MFI-cancer-meta-MR-DML`）

```
MFI-cancer-meta-MR-DML/
├── README.md                 # 项目英文主页：介绍、目录、数据来源、复现入口、DOI 徽章
├── LICENSE                   # MIT 开源协议（开放、宽松，最适合科研代码）
├── .gitignore                # 屏蔽规则：大包、缓存、密钥、本地环境一律不进库
├── R_analysis/               # 全部 R 代码
│   ├── 01_meta_analysis.R            # meta 分析：合并效应、异质性、留一法、偏倚、亚组 RVE
│   ├── 02_forward_MR.R               # 正向 2SMR（MFI/AS → 癌症）
│   ├── 03_reverse_MR.R               # 反向 2SMR（癌症 → MFI），已去除硬编码 token
│   ├── 04_sensitivity_FDR_MDOR.R     # 敏感性、BH-FDR 校正、MDOR 计算
│   ├── 05_forest_plot.R              # 森林图/漏斗图
│   └── renv.lock                      # R 环境锁（真实生成，R 4.5.1，含全部依赖版本）
├── Python_DML/               # 全部 Python 代码
│   ├── 01_dml_main.py                # DML 主分析（LassoCV + 随机森林，5 折交叉拟合）
│   ├── 02_figure3_plot.py            # 图 3 三联图与敏感性森林图
│   ├── requirements.txt              # pip 包版本锁定
│   └── environment.yml               # conda 环境锁定
├── output_results/           # 小体积结果表（按分析分目录）
│   ├── meta/  .gitkeep
│   ├── MR/    .gitkeep
│   └── DML/   .gitkeep
├── docs/
│   ├── run_instruction.md            # 英文复现说明（R 与 Python 两部分）
│   └── 代码归档与DOI操作指南.md       # 本文件
└── manuscript_supplement/    # 补充表 CSV（.gitkeep 占位，后续放入）
```

每个文件/文件夹的作用：

| 文件/文件夹 | 作用 |
| --- | --- |
| `README.md` | 别人打开仓库看到的第一页，说明项目是什么、怎么复现、数据在哪、如何引用 |
| `LICENSE` | 法律层面声明别人可以怎么用你的代码；MIT 允许自由使用、修改、引用，只需保留版权声明 |
| `.gitignore` | 告诉 Git 哪些文件"永远不要上传"，是防止泄密和仓库膨胀的关键 |
| `R_analysis/` | 按数字前缀编号，保证运行顺序清晰；`renv.lock` 锁定 R 包版本 |
| `Python_DML/` | DML 分析与绘图；`requirements.txt`/`environment.yml` 锁定 Python 包版本 |
| `output_results/` | 只放**小体积**结果表（CSV）；大图和原始数据不进库 |
| `docs/` | 复现说明与操作指南 |
| `manuscript_supplement/` | 论文补充表的 CSV 版本，便于审稿人核对 |
| `.gitkeep` | Git 不跟踪空文件夹，用这个占位文件让空目录也能被提交 |

### 1.2 `.gitignore` 已包含的屏蔽规则（要点）

- **原始 GWAS 大包**：`*.gz`、`*.zip`、`*.tar`、`*sumstats*` 等，全部屏蔽。
- **密钥/令牌**：`.env`、`*token*`、`*secret*`、`*credential*`，防止 OpenGWAS token 泄露。
- **R 缓存**：`.RData`、`.Rhistory`、`.Rproj.user/`、`*.Rout`，以及 `renv/library/` 本地库目录（保留 `renv.lock`，但不提交本地安装的包）。
- **Python 缓存与环境**：`__pycache__/`、`.venv/`、`venv/`、`.conda/`、`.ipynb_checkpoints/`、`*.pkl`、`*.joblib`。
- **系统/编辑器临时文件**：`.DS_Store`、`Thumbs.db`、`~$*`、`.idea/`、`.vscode/`。
- **输出目录里的大图**：`output_results/**/*.png/.tiff/.pdf` 屏蔽（保留小表格）。

### 1.3 R 环境锁定（renv）完整操作

本仓库的 `R_analysis/renv.lock` 已用你电脑上的 R 4.5.1 真实生成。以下是完整命令，供日后更新或在新机器上理解流程。

**在 RStudio 控制台里运行（首次，恢复环境）：**

```r
install.packages("renv")                       # 若未安装 renv
setwd("D:/研究生科研小组/小论文/MFI-cancer-meta-MR-DML")  # 改成你的仓库路径
renv::activate()                               # 激活项目环境
renv::restore()                                # 按 renv.lock 安装全部锁定版本的包
```

**生成/更新锁文件（以后加了新包才需要）：**

```r
renv::init()        # 首次初始化，自动扫描并快照已用的包
renv::snapshot()    # 增删包之后，更新 renv.lock
renv::status()      # 检查当前库与锁文件是否一致
```

> 本仓库锁文件已记录的关键包：`meta 8.5-0`、`metafor 5.0-1`、`robumeta 2.1`、`clubSandwich 0.7.0`、`readxl 1.5.0`、`data.table 1.17.8`、`remotes 2.5.0`，以及 GitHub 上的 `TwoSampleMR`、`ieugwasr`。

### 1.4 Python 环境导出（conda + pip）完整命令

**方式一：conda（推荐，已提供 `environment.yml`）**

```bat
:: 在 Anaconda Prompt 中，进入 Python_DML 目录
cd /d "D:\研究生科研小组\小论文\MFI-cancer-meta-MR-DML\Python_DML"

:: 按锁文件一键创建环境
conda env create -f environment.yml
conda activate dml-repro

:: 以后导出/更新环境的命令
conda env export --no-builds > environment.yml
```

**方式二：pip（已提供 `requirements.txt`）**

```bat
cd /d "D:\研究生科研小组\小论文\MFI-cancer-meta-MR-DML\Python_DML"
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt

:: 以后导出/更新
pip freeze > requirements.txt
```

> 锁定版本（Python 3.12）：numpy 2.1.2、pandas 2.2.3、scipy 1.14.1、scikit-learn 1.9.0、matplotlib 3.11.1、seaborn 0.13.2。

### 1.5 模板文件位置

- 英文 `README.md`、`LICENSE`、`.gitignore`：仓库根目录，已写好。
- 英文 `docs/run_instruction.md`：已写好，分 R 与 Python 两部分，含逐步命令和常见报错排查。
- 提交前只需把 `README.md`、`LICENSE` 里的占位信息（`[Author Name]`、`[Full manuscript title]`、Zenodo 编号）替换成真实内容。

---

## 模块 2　上传 GitHub：两种方案

### 为什么用 GitHub

GitHub 是**在线开发仓库**，负责版本管理、随时更新、多人协作，并提供一个可访问的代码链接。它解决"代码放哪、怎么追踪修改"的问题，但**本身不保证永久不变、也不自动给可引用 DOI**，所以正式归档还要靠模块 3 的 Zenodo。

### 方案 A（新手优先）：GitHub Desktop 图形界面

**第 1 步：注册账号并安装软件**
1. 浏览器打开 <https://github.com>，用邮箱注册一个账号（免费），验证邮箱。
2. 打开 <https://desktop.github.com>，下载安装 **GitHub Desktop**。
3. 打开 GitHub Desktop，点 *Sign in to GitHub.com*，用刚注册的账号登录。

**第 2 步：在 GitHub 网站新建仓库**
1. 登录 GitHub 网页，右上角点 **+ → New repository**。
2. *Repository name* 填规范英文名，例如 `MFI-cancer-meta-MR-DML`。
3. 选 **Private**（私有，投稿阶段不公开）或 **Public**（公开）。
4. 勾选 **Add a README file** 可留可不选（本地已有 README），*.gitignore* 和 license 本地已有，可都不选。
5. 点 **Create repository**。

**第 3 步：把仓库克隆到本地**
1. 回到 GitHub Desktop，菜单 *File → Clone repository*。
2. 在列表里选中刚建的仓库，*Local path* 选择一个**空文件夹或新位置**，点 *Clone*。
3. 建议先克隆到一个临时位置，再把你 `MFI-cancer-meta-MR-DML` 里的**全部内容**（含 `.gitignore` 等隐藏文件）粘贴进这个克隆文件夹。

**第 4 步：粘贴文件并核对没有违禁内容**
1. 打开克隆文件夹，把仓库内容粘贴进去。
2. 确认看不到任何 `.gz` 原始数据、`.RData`、`.venv`、token 文件——`.gitignore` 会自动拦截它们。
3. 回到 GitHub Desktop，左侧 *Changes* 列表里应只出现代码、CSV、文档，**仔细检查没有 token、没有大包**。

**第 5 步：提交（Commit）**
1. 在左下角 *Summary* 框填写提交说明，例如 `Initial reproducible code archive`。
2. 点 **Commit to main**。提交只是把改动记录在本地。

**第 6 步：推送（Push）到 GitHub**
1. 点顶部 **Push origin**（或 *Repository → Push*）。
2. 刷新 GitHub 网页，就能看到全部文件。以后每次改完代码，重复"写 Summary → Commit → Push"。

### 方案 B（备选）：Git 命令行

先安装 Git for Windows（<https://git-scm.com/download/win>），在仓库文件夹里右键选 *Git Bash Here*，然后：

```bash
# 1. 初始化本地仓库（在 MFI-cancer-meta-MR-DML 文件夹内）
git init
git branch -M main

# 2. 配置身份（仅首次）
git config --global user.name "Your Name"
git config --global user.email "you@example.com"

# 3. 查看将被提交的文件，确认没有大包/密钥（重点检查！）
git status

# 4. 暂存并提交
git add .
git commit -m "Initial reproducible code archive"

# 5. 关联 GitHub 远程仓库（URL 换成你在网站建好的仓库地址）
git remote add origin https://github.com/你的用户名/MFI-cancer-meta-MR-DML.git

# 6. 推送
git push -u origin main
```

> 以后日常更新：`git add .` → `git commit -m "说明"` → `git push`。

---

## 模块 3　对接 Zenodo 获取永久 DOI

### 3.1 为什么必须有 Zenodo DOI，GitHub 链接为何不能替代

- **GitHub 仓库会变**：你会持续 commit、甚至删除或设为私有，链接内容不固定，审稿人和读者无法保证多年后看到的还是"投稿那一版"。
- **Zenodo 是永久归档库**：由欧洲核子研究组织（CERN）运营，专门做科研成果长期保存。它把某一时刻的仓库**冻结成不可更改的快照**，并分配一个**永久、可引用的 DOI**。
- 标准做法是：**GitHub 负责日常开发，发布一个 Release 后由 Zenodo 抓取该 Release 生成 DOI**。论文里引用的是 Zenodo DOI，而不是会变动的 GitHub 地址。

### 3.2 分步操作

**第 1 步：注册/登录 Zenodo**
1. 打开 <https://zenodo.org>，右上角 *Sign up*，可直接用 GitHub 账号或机构邮箱登录（推荐用 GitHub 登录，方便授权）。

**第 2 步：授权 Zenodo 访问 GitHub**
1. 登录后打开 <https://zenodo.org/account/settings/github/>。
2. 点 **Connect GitHub**，按提示授权 Zenodo 读取你的 GitHub 账号。

**第 3 步：开启仓库同步开关**
1. 授权后，同一页面会列出你的 GitHub 仓库。
2. 找到 `MFI-cancer-meta-MR-DML`，把它右侧的开关**拨到 ON**。

**第 4 步：在 GitHub 发布一个 Release（这是归档的对象）**
1. 回到 GitHub 仓库网页，右侧点 **Releases → Create a new release**。
2. *Choose a tag* 输入版本号，建议 `v1.0.0`（语义化版本，投稿版用 1.0.0）。
3. *Release title* 填 `Reproducible code for manuscript submission`。
4. 可在描述里写明对应稿件，点 **Publish release**。

**第 5 步：Zenodo 自动生成 DOI**
1. 发布 Release 后等 1–2 分钟，Zenodo 会自动抓取并生成一条记录。
2. 回到 Zenodo，在 *Upload* 里能看到该记录，页面上有形如 `10.5281/zenodo.XXXXXXX` 的 **DOI**。
3. 点进记录，补全标题、作者、摘要、关键词、关联论文等元数据（这些让归档更规范、可被检索）。

**第 6 步：把 DOI 徽章添加到 README**
1. 在 Zenodo 记录页右侧 *DOI badge* 处复制 Markdown 代码。
2. 替换 `README.md` 顶部的占位徽章，把 `10.5281/zenodo.XXXXXXX` 改成真实编号。
3. 提交并推送，再发布一个小版本也可；DOI 已永久有效，不会因后续更新而失效。

> 关于"版本与 DOI"：每个新 Release 会生成一个**版本 DOI**，同时 Zenodo 会给整个项目一个**概念 DOI（Concept DOI）**，它始终指向最新版本。论文里引用具体版本 DOI 更严谨；若希望读者总是看到最新版，可引用概念 DOI。

---

## 模块 4　论文 Data and Code Availability 写法

### 版本一：返修阶段（尚未接收、还没有正式 DOI）

> The R code for the meta-analysis and bidirectional two-sample Mendelian randomization analyses, and the Python code for the summary-level double machine learning analyses, are available in a public GitHub repository at [GitHub repository URL]. The repository includes environment lock files (`renv.lock`, `environment.yml`, and `requirements.txt`) and step-by-step reproduction instructions. A permanent, citable archive with a DOI will be deposited in Zenodo upon acceptance, and the DOI will be provided in the final version. All GWAS summary statistics were obtained through the OpenGWAS API and are not redistributed here because of their data-use terms; the scripts retrieve the required data automatically.

### 版本二：接收校样阶段（已有 DOI，正式版）

> The complete analysis code supporting this study is publicly available and permanently archived in Zenodo under the DOI **10.5281/zenodo.XXXXXXX** (replace with the real DOI), and can also be browsed on GitHub at [GitHub repository URL]. The repository contains the R scripts for the meta-analysis and bidirectional two-sample Mendelian randomization, the Python scripts for the summary-level double machine learning, environment lock files (`renv.lock`, `environment.yml`, and `requirements.txt`), and step-by-step reproduction instructions. GWAS summary statistics were obtained from the OpenGWAS API and are not redistributed because of their data-use terms; all data can be retrieved automatically by the provided scripts. The observational meta-analysis input table is documented within the repository.

> 中文版（若期刊/学位论文要求中文）：本研究的全部分析代码已通过 Zenodo 永久归档并分配可引用 DOI（10.5281/zenodo.XXXXXXX），同时可在 GitHub 仓库（[URL]）浏览。仓库包含 meta 分析与双向两样本孟德尔随机化的 R 脚本、汇总水平双机器学习的 Python 脚本、环境锁定文件及逐步复现说明。GWAS 汇总数据通过 OpenGWAS 接口获取，受数据使用条款限制不在仓库中重新分发，脚本可自动调取所需数据。

---

## 附：避坑清单（提交前逐条核对）

1. **绝不上传原始 GWAS 大包**：`.gz` 等 sumstats 体积可达数百 MB 到数 GB，GitHub 单文件上限 100 MB，且多有数据许可限制。已由 `.gitignore` 屏蔽，代码通过 API 在线取数。
2. **绝不硬编码或提交 OpenGWAS token / API key / 密码**：仓库里的 `03_reverse_MR.R` 已改为从环境变量 `OPENGWAS_JWT` 读取；提交前用 `git status` 和搜索确认没有 `eyJ`、`token`、`.env` 等字样。一旦 token 曾被提交，应立即到 OpenGWAS 作废旧 token 并重新生成。
3. **不上传缓存与本地环境**：`.RData`、`.Rhistory`、`.Rproj.user/`、`renv/library/`、`__pycache__/`、`.venv/`、`.ipynb_checkpoints/` 都不进库；只提交锁文件，不提交整个包库。
4. **不上传模型/临时大文件**：`*.pkl`、`*.joblib`、`*.h5`、`*.tmp`、`~$*`、`Thumbs.db` 等。
5. **文件名与目录用规范英文、不用空格和特殊字符**，脚本用数字前缀保证运行顺序；中文路径只出现在你本机，仓库内引用尽量用相对路径。
6. **先私有、后公开**：投稿期间可把 GitHub 仓库设为 Private；需要审稿人查看时用 Zenodo 归档或生成 GitHub 临时访问链接，接收后再公开。
7. **Release 版本号要规范**：投稿归档用 `v1.0.0`，后续修改递增（v1.1.0、v2.0.0），不要反复覆盖同一 Release。
8. **DOI 引用要对应冻结版本**：论文写的是 Zenodo 版本 DOI，不要只写会变动的 GitHub 链接；README 里的 DOI 徽章占位要替换成真实编号。
9. **提交前做一次"干净机测试"**：按 `docs/run_instruction.md` 在新环境里从锁文件恢复并跑通，确认没有依赖缺失、没有写死的本机绝对路径导致的报错。
