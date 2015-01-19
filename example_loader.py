from models import db, Resource, User, Client
from indexer import index_search_elements
from fhir_parser import parse_resource
from fhir_spec import RESOURCES
import names
from vcf import VCFReader
from random import random
from functools import partial
from config import MAX_SEQ_PER_FILE
import os


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
    gender = 'female' if random() < 0.5 else 'male'
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
    print 'Created patient called %s'% full_name
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
    print 'Created procedure (sequencing lab)'
    return save_resource('Procedure', data)


def load_patients_by_samples(samples):
    return {sample: rand_patient() for sample in samples}

def load_lab_by_patients(patients):
    # patients is a key-value pair of sample and patient
    return {sample: rand_lab(patients[sample])
        for sample in patients.keys()}


def load_vcf_example(vcf_file):
    reader = VCFReader(filename=vcf_file)
    patients = load_patients_by_samples(reader.samples)
    db.session.commit()
    labs = load_lab_by_patients(patients)
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
            referenced_patient = patients[sample.sample]
            referenced_lab = labs[sample.sample]
            seq_data['patient'] = referenced_patient.get_reference()
            seq_data['source']['lab'] = referenced_lab.get_reference()
            variant = record.ID if record.ID is not None else 'anonymous variant'
            seq_data['text']['div']  = '<div>Genotype of %s is %s</div>'% (variant, reads)
            save_resource('Sequence', seq_data)
            print 'Created sequence at %s:%s-%s'% (record.CHROM, record.POS, record.end)
            count += 1

        if count > MAX_SEQ_PER_FILE:
            break


def init_superuser():
    superuser = User(email='super')
    db.session.add(superuser)
    global test_resource
    test_resource = partial(Resource, owner_id=superuser.email)  


def load_examples():
    from os import path
    basedir = path.dirname(path.abspath(__file__))
    init_superuser()
    for example_file in os.listdir(path.join(basedir, 'examples/vcf')):
        load_vcf_example(path.join(basedir, 'examples/vcf', example_file))
