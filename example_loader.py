from models import db, Resource, User, Client
from indexer import index_search_elements
from fhir_parser import parse_resource
from fhir_spec import RESOURCES
import names
from vcf import VCFReader
import random
from functools import partial
from config import MAX_SEQ_PER_FILE, CONDITION_TO_SEQ_RATIO
import json
import os


BASEDIR = os.path.dirname(os.path.abspath(__file__))

RELIABILITIES = ['questionable', 'ongoing', 'ok', 'calibrating', 'early']
CYP2C19_STAR_VARIANTS = [
  {
    'gene' : {
      'display': 'CYP2C19 Gene',
      'code': '124020',
      'system': 'http://www.ncbi.nlm.nih.gov/omim'
    }, 
    'value' : ['*1/*1', '*2/*17']
  }
  ]

SLCO1B1_SNPS = [
  {
    'id': 'rs4149056',
    'chrom': '12',
    'pos': 21178615,
    'value': ['T/T', 'C/C', 'C/T']
  }
]

CYP2C19_SNPS = [
  {'id': 'rs4244285', 'chrom': '10', 'pos': 94781859 },
  {'id': 'rs4986893', 'chrom': '10', 'pos': 94780653 },
  {'id': 'rs28399504', 'chrom': '10', 'pos': 94762706 },
  {'id': 'rs56337013', 'chrom': '10', 'pos': 94852738 },
  {'id': 'rs72552267', 'chrom': '10', 'pos': 94775453 },
  {'id': 'rs72558186', 'chrom': '10', 'pos': 94781999 },
  {'id': 'rs41291556', 'chrom': '10', 'pos': 94775416 },
  {'id': 'rs12248560', 'chrom': '10', 'pos': 94761900 },
]
  
INTERPRETATIONS = [
    {
        'code': 'L',
        'display': 'Below low normal',
        'system': 'http://hl7.org/fhir/vs/observation-interpretation'
    }, { 
        'code': 'IND',
        'display': 'Intermediate',
        'system': 'http://hl7.org/fhir/vs/observation-interpretation'
    }, { 
        'code': 'H',
        'display': 'Above high normal',
        'system': 'http://hl7.org/fhir/vs/observation-interpretation'
    }, { 
        'code': 'NEG',
        'display': 'Negative',
        'system': 'http://hl7.org/fhir/vs/observation-interpretation'
    }, { 
        'code': 'POS',
        'display': 'Positive',
        'system': 'http://hl7.org/fhir/vs/observation-interpretation'
    }
]


def save_resource(resource_type, resource_data):
    '''
    save a resource to database and index its elements by search params
    '''
    valid, search_elements = parse_resource(resource_type, resource_data)
    assert valid
    resource = test_resource(resource_type, resource_data) 
    index_search_elements(resource, search_elements)
    db.session.add(resource)
    return resource


def rand_patient():
    '''
    generate random resource and index its elements by search params
    '''
    gender = 'female' if random.random() < 0.5 else 'male'
    first_name = names.get_first_name(gender=gender)
    last_name = names.get_last_name()
    full_name = '%s %s'% (first_name, last_name)
    data = {
        'resourceType': 'Patient',
        'text': {
            'status': 'generated',
            'div': "<div><p>%s</p></div>"% full_name
        },
        'name': [{'text': full_name}],
        'gender': {
                'text': gender,
                'coding': [{
                    'code': 'F' if gender == 'female' else 'M',
                    'system': 'http://hl7.org/fhir/v3/AdministrativeGender'}]
        }
    }
    print 'Created Patient called %s'% full_name
    return save_resource('Patient', data)



def rand_lab(patient):
    data = {
        'resourceType': 'Procedure',
        'text': {
            'status': 'generated',
            'div': '<div>DNA Sequencing lab</div>'
        },
        'subject': patient.get_reference(),
        'type':  {
            'text': 'Sequencing lab',
            'coding': [{
                'code': 'SequencingLab',
                'system': 'http://genomics.smartplatforms.org/dictionary#sequencinglab'
            }]
        }
    }
    print 'Created Procedure (Sequencing Lab)'
    return save_resource('Procedure', data)

    
#def rand_rx(patient):
def clopdigorel_rx(patient, dose):
    data = {
      "resourceType": "MedicationPrescription",
      "text": {
        "status": "generated",
        "div": "<div>clopidogrel " + dose + "mg, 1/day</div>"
      },
      "status": "active",
      "patient": patient.get_reference(),
      "medication": {
        "reference": "Medication/clopidogrel"
      },
      "dosageInstruction": [
          {
            "timingSchedule": {
              "repeat": {
                "frequency": 1,
                "duration": 1,
                "units": "d"
              }
            },
            "doseQuantity": {
              "value": float(dose),
              "units": "mg",
              "system": "http://unitsofmeasure.org",
              "code": "mg"
            }
          }
        ]
    }

    print 'Created Medication Prescription'
    return save_resource('MedicationPrescription', data)

def simvastatin_rx(patient, dose):
    data = {
      "resourceType": "MedicationPrescription",
      "text": {
        "status": "generated",
        "div": "<div>simvastatin " + dose + "mg, 1/day</div>"
      },
      "status": "active",
      "patient": patient.get_reference(),
      "medication": {
        "reference": "Medication/simvastatin"
      },
      "dosageInstruction": [
          {
            "timingSchedule": {
              "repeat": {
                "frequency": 1,
                "duration": 1,
                "units": "d"
              }
            },
            "doseQuantity": {
              "value": float(dose),
              "units": "mg",
              "system": "http://unitsofmeasure.org",
              "code": "mg"
            }
          }
        ]
    }

    print 'Created Medication Prescription'
    return save_resource('MedicationPrescription', data)
    
def load_patients_by_samples(samples):
    return {sample: rand_patient() for sample in samples}


def load_labs_by_patients(patients):
    # patients is a key-value pair of sample and patient
    return {sample: rand_lab(patients[sample])
        for sample in patients.keys()}

#def load_meds_by_patients(patients):
#    return {sample: rand_rx(patients[sample])
#        for sample in patients.keys()}
        
'''def load_star_variants_by_patients(patients, star_variant_list):
    return {sample: make_star_variant(patients[sample], star_variant_list)
        for sample in patients.keys()}'''

def rand_conditions(patient):
    '''
    randomly assign a set of conditions to a poor patient
    '''
    conditions = random.sample(available_conditions,  
                            random.randint(0, len(available_conditions)))
    ret = []
    for cond_tmpl in conditions:
        print cond_tmpl['code']['text']
        condition = dict(cond_tmpl)
        condition['subject'] = patient.get_reference()
        ret.append(save_resource('Condition', condition))
        print 'Created condition %s'% condition['code'].get('text', '')

    return ret

'''def make_star_variant(patient, star_variant_list):
    random_result = random.choice(star_variant_list);
    observation = {
        'resourceType': 'Observation',
        'subject': patient.get_reference(),
        'name': {
            'coding': [random_result['gene']]
        },
        'valueString': random.choice(random_result['value']),
        'status': 'final',
        'reliability': random.choice(RELIABILITIES)
    }
    print 'Created Observation (Genetic Star Variant Observation)'
    print observation
    return save_resource('Observation', observation)'''
    
def make_star_variant(patient, gene, value):
    observation = {
        'resourceType': 'Observation',
        'subject': patient.get_reference(),
        'name': {
            'coding': [gene]
        },
        'valueString': value,
        'status': 'final',
        'reliability': random.choice(RELIABILITIES)
    }
    print 'Created Observation (Genetic Star Variant Observation)'
    print observation
    return save_resource('Observation', observation)

    
def make_observation(condition, sequence, patient, seq_id):
    observation = {
        'resourceType': 'Observation',
        'extension': [
            {
                'url': 'http://genomics.smartplatforms.org/dictionary/GeneticObservation#assessedCondition',
                'valueReference': condition.get_reference()
            }, {
                'url': 'http://genomics.smartplatforms.org/dictionary/GeneticObservation#variantGenotype',
                'valueReference': sequence.get_reference()
            }
        ],
        'subject': patient.get_reference(),
        'name': {
            'coding': [{
                'display': 'Genetic Observation',
                'code': 'GeneticObservation',
                'system': 'http://genomics.smartplatforms.org/dictionary'
            }]
        },
        'interpretation': random.choice(INTERPRETATIONS),
        'status': 'final',
        'reliability': random.choice(RELIABILITIES)
    }
    if seq_id is not None:
        observation['extension'].append({
                'url': 'http://genomics.smartplatforms.org/dictionary/GeneticObservation#variantId',
                'valueString': seq_id
        })

    print 'Created Observation (Genetic Observation)'
    return save_resource('Observation', observation)



def load_conditions_by_patients(patients):
    return {sample: rand_conditions(patients[sample])
        for sample in patients.keys()}

        
'''def create_snp_result(patient, lab, snp_list):
    snp = random.choice(snp_list)
    sequence_tmpl = {
        'text': {'status': 'generated'},
        'resourceType': 'Sequence',
        'type': 'dna',
        'chromosome': snp['chrom'],
        'startPosition': snp['pos'],
        'endPosition': snp['pos'],
        'assembly': 'GRCh37',
        'source': {'sample': 'somatic'},
        'snp': snp['id']
    }
    
    seq_data = dict(sequence_tmpl)
    reads = random.choice(snp['value'])
    seq_data['read'] = reads.split('/')
    # links sequence to patient and lab
    referenced_patient = patient
    referenced_lab = lab
    seq_data['quality'] = 68.33
    seq_data['patient'] = patient.get_reference()
    seq_data['source']['lab'] = referenced_lab.get_reference()
    # get name of the variant
    variant_id = snp['id']
    variant = variant_id if variant_id is not None else 'anonymous variant'
    seq_data['text']['div']  = '<div>Genotype of %s is %s</div>'% (variant, reads)
    print seq_data
    sequence = save_resource('Sequence', seq_data)
    print 'Created SNP at %s:%s  %s'% (snp['chrom'], snp['pos'], seq_data['text']['div'])'''
    
def create_snp_result(patient, lab, snp, reads):
    sequence_tmpl = {
        'text': {'status': 'generated'},
        'resourceType': 'Sequence',
        'type': 'dna',
        'chromosome': snp['chrom'],
        'startPosition': snp['pos'],
        'endPosition': snp['pos'],
        'assembly': 'GRCh37',
        'source': {'sample': 'somatic'},
        'snp': snp['id']
    }
    
    seq_data = dict(sequence_tmpl)
    seq_data['read'] = reads.split('/')
    # links sequence to patient and lab
    referenced_patient = patient
    referenced_lab = lab
    seq_data['quality'] = 68.33
    seq_data['patient'] = patient.get_reference()
    seq_data['source']['lab'] = referenced_lab.get_reference()
    # get name of the variant
    variant_id = snp['id']
    variant = variant_id if variant_id is not None else 'anonymous variant'
    seq_data['text']['div']  = '<div>Genotype of %s is %s</div>'% (variant, reads)
    print seq_data
    sequence = save_resource('Sequence', seq_data)
    print 'Created SNP at %s:%s  %s'% (snp['chrom'], snp['pos'], seq_data['text']['div'])
    
'''def load_snps_by_patients(patients, labs, snp_list):
    return {sample: create_snp_result(patients[sample], labs[sample], snp_list)
        for sample in patients.keys()}'''

def create_patient(is_clop_star, is_clop_snp, is_slco1b1_normal, is_slco1b1_on_med):
    patient = rand_patient()
    db.session.commit()
    if is_clop_star or is_clop_snp:
        clopdigorel_rx(patient, '75')
        if is_clop_star:
            make_star_variant(patient, CYP2C19_STAR_VARIANTS[0]['gene'], '*1/*1')
        else:
            lab = rand_lab(patient)
            db.session.commit()
            # Profile for a *2/*3 patient
            create_snp_result(patient, lab, CYP2C19_SNPS[0], 'A/G')
            create_snp_result(patient, lab, CYP2C19_SNPS[1], 'A/G')
            create_snp_result(patient, lab, CYP2C19_SNPS[2], 'A/A')
            create_snp_result(patient, lab, CYP2C19_SNPS[3], 'C/C')
            create_snp_result(patient, lab, CYP2C19_SNPS[4], 'G/G')
            create_snp_result(patient, lab, CYP2C19_SNPS[5], 'T/T')
            create_snp_result(patient, lab, CYP2C19_SNPS[6], 'T/T')
            create_snp_result(patient, lab, CYP2C19_SNPS[7], 'C/C')
    else:
        lab = rand_lab(patient)
        db.session.commit()
        if is_slco1b1_on_med:
            simvastatin_rx(patient, '20')
            create_snp_result(patient, lab, SLCO1B1_SNPS[0], 'C/C')
        else:
            create_snp_result(patient, lab, SLCO1B1_SNPS[0], 'T/T')
        
def load_vcf_example(vcf_file):
    reader = VCFReader(filename=vcf_file)
    patients = load_patients_by_samples(reader.samples)
    db.session.commit()
    meds = load_meds_by_patients(patients)
    load_star_variants_by_patients(patients, CYP2C19_STAR_VARIANTS)
    labs = load_labs_by_patients(patients)
    db.session.commit()
    conditions = load_conditions_by_patients(patients)
    db.session.commit()
    load_snps_by_patients(patients, labs, SLCO1B1_SNPS)
    db.session.commit()
    count = 0
    for record in reader:
        sequence_tmpl = {
            'text': {'status': 'generated'},
            'resourceType': 'Sequence',
            'type': 'dna',
            'chromosome': record.CHROM,
            'startPosition': record.POS,
            'endPosition': record.end,
            'assembly': 'GRCh37',
            'source': {'sample': 'somatic'}
        }
        for sample in record.samples:
            sample_id = sample.sample
            reads = sample.gt_bases
            if '/' in reads:
                delimiter = '/'
            elif '|' in reads:
                delimiter = '|'
            else:
                delimiter = '.'
            seq_data = dict(sequence_tmpl)
            seq_data['read'] = reads.split(delimiter)
            # get genotype quality 
            if 'GQ' in dir(sample.data):
                seq_data['quality'] = sample.data.GQ
            # links sequence to patient and lab
            referenced_patient = patients[sample_id]
            referenced_lab = labs[sample_id]
            seq_data['patient'] = referenced_patient.get_reference()
            seq_data['source']['lab'] = referenced_lab.get_reference()
            # get name of the variant
            variant_id = record.ID
            seq_data['snp'] = variant_id
            variant = variant_id if variant_id is not None else 'anonymous variant'
            seq_data['text']['div']  = '<div>Genotype of %s is %s</div>'% (variant, reads)
            #print seq_data
            sequence = save_resource('Sequence', seq_data)
            print 'Created Sequence at %s:%s-%s  %s'% (record.CHROM, record.POS, record.end, seq_data['text']['div'])

            # randomly link a DNA sequence to conditions that the user has
            if (len(conditions[sample_id]) > 0 and
                random.random() <= CONDITION_TO_SEQ_RATIO):
                make_observation(random.choice(conditions[sample_id]),
                                sequence,
                                referenced_patient,
                                variant_id)
            count += 1

        if count >= MAX_SEQ_PER_FILE:
            break


def load_condition_from_file(path):
    abspath = os.path.join(BASEDIR, 'examples/conditions', path)
    with open(abspath) as condition_f:
        return json.loads(condition_f.read())


def init_conditions():
    condition_dir = os.path.join(BASEDIR, 'examples/conditions')
    global available_conditions
    available_conditions = map(load_condition_from_file, os.listdir(condition_dir))


def init_superuser():
    superuser = User(email='super')
    db.session.add(superuser)
    global test_resource
    test_resource = partial(Resource, owner_id=superuser.email)  


def load_examples():
    init_superuser()
    init_conditions()
    create_patient(True, False, False, False)
    create_patient(False, True, False, False)
    create_patient(False, False, True, False)
    create_patient(False, False, False, True)
    #for example_file in os.listdir(os.path.join(BASEDIR, 'examples/vcf')):
    #    load_vcf_example(os.path.join(BASEDIR, 'examples/vcf', example_file))
