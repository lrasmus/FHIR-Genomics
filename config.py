APP_CONFIG = {
	'SQLALCHEMY_DATABASE_URI': "sqlite:///fhir.db"
}

FHIR_SPEC_DIR = ''
# the example resource parser right now is still pretty slow
# also, due to the fact that SQLAlchemy is not very good at bulk insert,
# if the vcf file is too large, the example loader is gonna run for decades
MAX_SEQ_PER_FILE = 50
# ratio between Conditions and Sequences, the example loader uses this to
# randomly create GeneticObservation, which associates the two resources.
CONDITION_TO_SEQ_RATIO = 0.1
