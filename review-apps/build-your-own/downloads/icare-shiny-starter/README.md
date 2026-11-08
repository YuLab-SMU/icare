# ICARE Shiny starter kit

Guide: https://leonardfei.github.io/Icare_test/review-apps/build-your-own/

This kit contains no fitted models, patient data or account credentials. Add your own trained objects before running it.
Reference version: icare 1.1.0, YuLab-SMU/icare@7ff676a556721affeb6da6c36b3bb60fc24ca468; tested with R 4.5.1.
Use the R environment in which you trained your models. Model objects and learner packages must be compatible with this version.

1. Choose apps/classification, apps/subtyping or apps/survival.
2. Add your own model.rds: a Train_Model, Subtyping or PrognosiX object, respectively.
   The classification template reads the caret training schema from train.models and also requires preprocess.rds saved during training.
   If no external preprocessing was used, saveRDS(NULL, "apps/classification/preprocess.rds").
3. From the kit root, run Rscript make_example.R survival (replace the application type as needed). Review the generated artificial profiles and their units.
4. Edit the title, study description, variable descriptions, units and citation in app.R. Declare any additional learner packages with library().
5. Run Rscript check_app.R survival, then manually test shiny::runApp("apps/survival").
   Subtyping checks nmf by default; use a third argument of lpa or kmeans to test another fitted method.
   Classification checks the first available model. Test every additional method you intend to expose.
6. Register your own shinyapps.io account in your local R console. Never include tokens or secrets in shared files.
7. Run Rscript deploy.R survival YOUR_ACCOUNT my-survival-app.
8. Open the hosted application and test prediction, supported CSV upload and download. Inspect the server logs. A page loading successfully does not prove prediction works.

RDS objects may contain training records. Review their contents and data permissions before deployment. Do not put model files in www/ or a public GitHub repository.
www/example.csv contains artificial marginal-quantile profiles for interface testing. It is not a validation dataset and may contain implausible combinations.
The survival template includes ranger for the example surv.ranger model. Add the appropriate dependencies and revalidate when using a different learner.
Train models and fit preprocessing outside app.R. Do not refit scaling or imputation on each new input.
