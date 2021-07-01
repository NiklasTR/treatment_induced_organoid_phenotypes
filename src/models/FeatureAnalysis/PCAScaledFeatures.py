from OrganoidFeatures import OrganoidFeatures
import DrugEffects
import LineDifferences
import Utils
import Config
import CreateOrganoidCutouts

# Calculate drug-induced phenotypes
lines = Utils.get_all_lines("human")
print(lines)
#LineDifferences.process_dmso_organoids_no_filter("human")
LineDifferences.process_all_organoids("human")

