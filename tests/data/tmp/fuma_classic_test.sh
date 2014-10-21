
fuma \
	-a "hg18:test_FusionDetectionExperiment.TestFusionDetectionExperiment.bed" \
	\
	-s  "ExperimentA:chimerascan:test_FusionDetectionExperiment.TestFusionDetectionExperiment.bedpe" \
	    "ExperimentB:chimerascan:test_FusionDetectionExperiment.TestFusionDetectionExperiment.bedpe" \
	    "ExperimentC:chimerascan:test_FusionDetectionExperiment.TestFusionDetectionExperiment.bedpe" \
	\
	-l "ExperimentA:hg18" \
	   "ExperimentB:hg18" \
	   "ExperimentC:hg18" \
	\
	-o "fuma_classic_"

# ----------------------------------------------------------------------

