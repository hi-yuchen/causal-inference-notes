# 从 OLS 到因果识别

[在线阅读中文讲义](https://hi-yuchen.github.io/causal-inference-notes/) · [第一篇原文](01-probability.md)

概率论与统计基础、从 OLS 到因果识别、RCT、RDD、DID、IV，以及控制变量、固定效应和标准误专题。七个章节文件均可直接在 GitHub 阅读。离线阅读请用浏览器打开 `index.html`。

## 构建与练习

修改编号 Markdown 文件后，运行 `python3 build.py` 重新生成网页。需要 Python 3.9+ 与 Pandoc 3。GitHub Pages 发布 main 分支根目录。

- `python3 ames.py` 枚举 Horvitz–Thompson 案例的 190 种样本，复算抽样设计比较。
- `designs.do` 包含 Stata 16+ 的五类教学比较及系数断言；本版未在 Stata 环境运行。
