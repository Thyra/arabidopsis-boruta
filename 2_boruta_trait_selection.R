# This script performs a Boruta trait selection on the wide-shaped table in intermediary_data/1745AJ_Phenotyping_nafixed.csv
# Results are saved unfixed and roughFixed (with or without Tentative decisions)

library('Boruta')
set.seed(17022019)

data.nafixed = read.csv("data/intermediary/1745AJ_Phenotyping_nafixed.csv")
data.nafixed$Treatment = as.factor(data.nafixed$Treatment)

boruta_result = Boruta(Treatment ~ ., data=data.nafixed, maxRuns=1500, doTrace=2)
trait_importance = attStats(boruta_result)
trait_importance = cbind(trait = rownames(trait_importance), trait_importance)
# Some results from playing with this:
#  - order of most important traits varies a bit each time this is run but median trait importance stays quite stable (< 1% difference between runs)
write.csv(trait_importance, "data/results/trait_selection_result.csv", row.names=F)

trait_importance_roughfixed = attStats(TentativeRoughFix(boruta_result))
trait_importance_roughfixed = cbind(trait = rownames(trait_importance_roughfixed), trait_importance_roughfixed)
write.csv(trait_importance_roughfixed, "data/results/trait_selection_result_roughFixed.csv", row.names=F)
