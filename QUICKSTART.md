# 🚀 PharmaAgentBench 快速开始指南

## 一、环境准备（5 分钟）

### 1. 安装 Python 依赖

```bash
# 使用 pip 安装
pip install -r requirements.txt

# 或者使用 conda（推荐）
conda create -n pharma python=3.10
conda activate pharma
pip install -r requirements.txt
```

**注意**: RDKit 的安装可能需要额外步骤：

```bash
# 使用 conda 安装 RDKit (推荐)
conda install -c conda-forge rdkit

# 或使用 pip (可能在某些系统上有问题)
pip install rdkit
```

### 2. 配置 API 密钥

#### 方式一：编辑配置文件
编辑 `config/api_keys.yaml`，填入你的 API 密钥：

```yaml
deepseek:
  api_key: "sk-your-deepseek-key"  # 替换为真实密钥
  
qwen:
  api_key: "sk-your-dashscope-key"  # 替换为真实密钥
  
gemini:
  api_key: "your-gemini-key"        # 替换为真实密钥
```

#### 方式二：使用环境变量（更安全）
```bash
# Windows PowerShell
$env:DEEPSEEK_API_KEY="sk-your-key"
$env:DASHSCOPE_API_KEY="sk-your-key"
$env:GEMINI_API_KEY="your-key"

# Linux/Mac
export DEEPSEEK_API_KEY="sk-your-key"
export DASHSCOPE_API_KEY="sk-your-key"
export GEMINI_API_KEY="your-key"
```

### 3. 获取 API 密钥

| 平台 | 注册地址 | 说明 |
|------|---------|------|
| DeepSeek | https://platform.deepseek.com/ | 国产大模型，性价比高 |
| Qwen (DashScope) | https://dashscope.aliyun.com/ | 阿里云通义千问 |
| Gemini | https://ai.google.dev/ | Google AI |

---

## 二、运行演示（2 分钟）

### 快速体验
```bash
python main.py --demo
```

这将运行一个预设的分子生成任务，展示完整流程。

### 预期输出
```
    ╔═══════════════════════════════════════════════════════════════╗
    ║           🧬 PharmaAgentBench 医药大模型评测系统 🧬            ║
    ╚═══════════════════════════════════════════════════════════════╝

🔑 检查 API 密钥配置...
  ✓ DEEPSEEK: 已配置
  ✓ QWEN: 已配置
  ✓ GEMINI: 已配置

🚀 初始化 AgentScope 框架...
✓ AgentScope 初始化完成

📋 任务配置:
  类型: molecule_generation
  目标: 设计一个针对 EGFR 突变体的小分子抑制剂
  模型: deepseek, qwen, gemini

🤖 启动多智能体系统...
⚙️  执行评测任务...

[... 执行过程 ...]

📊 评测报告摘要
...
💾 报告已保存至: reports/report_molecule_generation_20251210_143022.html
```

---

## 三、交互式使用（5 分钟）

### 启动交互模式
```bash
python main.py
```

### 选择任务类型

```
请选择任务类型:
  1. 分子生成 (molecule_generation)
  2. 文献信息提取 (literature_extraction)
  3. ADMET 预测 (admet_prediction)
  4. 退出

请输入选项 (1-4): 1
```

### 输入任务描述

```
请输入任务目标描述: 设计一个 CDK4/6 抑制剂，用于治疗乳腺癌
```

### 等待结果

系统将：
1. ✅ 调用 3 个大模型生成分子
2. ✅ 使用 RDKit 验证化学合法性
3. ✅ 检查类药性和毒性
4. ✅ 生成 HTML 可视化报告

---

## 四、查看报告（3 分钟）

### 报告位置
```
reports/
└── report_molecule_generation_20251210_143022.html
```

### 在浏览器中打开
双击 HTML 文件，或在命令行中：

```bash
# Windows
start reports\report_*.html

# Mac
open reports/report_*.html

# Linux
xdg-open reports/report_*.html
```

### 报告内容包括：

1. **任务信息概览**
   - 任务类型、目标描述、评测时间
   
2. **模型对比表格**
   - 各模型的执行状态
   - 验证结果（通过/失败）
   - 警告和错误数量

3. **详细结果**
   - 每个模型生成的分子 SMILES
   - 分子性质（分子量、LogP、QED 等）
   - 设计理由和合成路线
   - 化学验证结果
   - 毒性警告

---

## 五、常见任务示例

### 1. 分子生成任务

```python
# 在交互模式下输入:
任务类型: 1 (molecule_generation)
目标描述: 设计一个针对 PD-L1 的小分子抑制剂，要求分子量 < 450
```

### 2. 文献提取任务

```python
任务类型: 2 (literature_extraction)
目标描述: 提取 BTK 抑制剂的临床试验数据和分子结构
```

### 3. ADMET 预测任务

```python
任务类型: 3 (admet_prediction)
目标描述: 预测 CC(C)Cc1ccc(cc1)C(C)C(O)=O 的 ADMET 性质
```

---

## 六、高级配置

### 自定义约束条件

编辑 `main.py` 中的 `demo_task`：

```python
demo_task = {
    "task_type": "molecule_generation",
    "target": "你的目标描述",
    "models": ["deepseek", "qwen", "gemini"],  # 可选择部分模型
    "constraints": {
        "max_mw": 450,          # 最大分子量
        "min_qed": 0.6,         # 最小 QED 评分
        "check_toxicity": True  # 是否检查毒性
    },
    "output_format": "html"
}
```

### 自定义 Prompt

编辑 `prompts/molecule_generation.yaml`：

```yaml
system_prompt: |
  你是一个专业的药物化学AI助手...
  
  # 在这里添加你的自定义指令
```

---

## 七、故障排除

### 问题 1: RDKit 安装失败

**解决方案**: 使用 conda 安装
```bash
conda install -c conda-forge rdkit
```

### 问题 2: API 调用失败

**检查清单**:
- ✅ API 密钥是否正确填写
- ✅ 网络是否可以访问对应平台
- ✅ API 额度是否用尽

**调试命令**:
```bash
python -c "from config import get_config; print(get_config().validate_api_keys())"
```

### 问题 3: AgentScope 初始化失败

**解决方案**: 检查版本
```bash
pip show agentscope
# 确保版本 >= 1.0.0
```

### 问题 4: 报告未生成

**检查**:
- ✅ `reports/` 目录是否存在
- ✅ 是否有写入权限

**手动创建目录**:
```bash
mkdir reports
```

---

## 八、下一步

### 扩展功能
- 🔧 添加更多评测指标（相似度分析、分子多样性）
- 🔧 集成更多大模型（Claude、GPT-4 等）
- 🔧 支持批量评测和基准测试
- 🔧 添加分子可视化（2D/3D 结构图）

### 学习资源
- [AgentScope 文档](https://doc.agentscope.io/)
- [RDKit 教程](https://www.rdkit.org/docs/)
- [DeepSeek API 文档](https://platform.deepseek.com/docs)

---

## 🎉 开始你的评测之旅！

```bash
python main.py --demo
```

有问题？请查看 [README.md](README.md) 或提交 Issue。
