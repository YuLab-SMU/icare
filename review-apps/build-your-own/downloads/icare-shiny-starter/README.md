# ICARE Shiny starter kit

教程：https://leonardfei.github.io/Icare_test/review-apps/build-your-own/

此模板不含模型、患者数据或账户凭据，不能在未添加自己的模型时直接预测。
基于 icare 1.1.0，YuLab-SMU/icare@7ff676a556721affeb6da6c36b3bb60fc24ca468；验证环境为 R 4.5.1。
在训练这些模型的 R 环境中操作。模型对象和算法后端须与该版本兼容。

1. 选择 apps/classification、apps/subtyping 或 apps/survival。
2. 放入自己的 model.rds，分别要求 Train_Model、Subtyping、PrognosiX 对象。
   classification 模板使用 caret 的 train.models 训练字段，还需训练阶段保存的 preprocess.rds。
   未做外部预处理时保存 saveRDS(NULL, "apps/classification/preprocess.rds")。
3. 根目录运行 Rscript make_example.R survival（可替换类型）。检查生成的人工样例单位与取值。
4. 编辑 app.R 的标题、研究描述、变量说明和单位；依实际算法添加 library() 声明。
5. 根目录运行 Rscript check_app.R survival，再用 shiny::runApp("apps/survival") 手工复核。
   分型默认检查 nmf，可用第三个参数指定 lpa 或 kmeans；只检查自己确实拟合的方法。
   分类自动检查首个可选模型；公开其他算法前也要逐一测试。
6. 在本机 R 控制台绑定自己的 shinyapps.io 账户；不要把 token/secret 放进任何共享文件。
7. 根目录运行 Rscript deploy.R survival YOUR_ACCOUNT my-survival-app。
8. 发布后重新用浏览器执行预测、上传及下载；查看云端日志。页面可打开不代表预测可用。

模型 RDS 可能包含原始训练记录。先审查数据和使用权限，不要放进 www/ 或公开 GitHub 仓库。
www/example.csv 是人工分位数组合，供界面测试，不是验证集，也不保证医学上合理。
生存模板的 ranger 后端适用于示例 surv.ranger 模型；使用其他 learner 时补齐其对应依赖并重新验证。
训练与预处理应在 app.R 之外完成。不要对每个新输入重新拟合标准化或插补步骤。
