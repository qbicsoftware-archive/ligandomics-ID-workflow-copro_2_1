from CTDopts.CTDopts import _InFile, CTDModel, args_from_file
import sys
import os
import subprocess
import re
import pandas as pd
from itertools import combinations



#### QBIC QPortal parameter transfer #####
pattern = re.compile('Q\w{4}[0-9]{3}[a-zA-Z]\w')

wf_dir = sys.argv[1]
ctd_params = args_from_file(wf_dir + '/WORKFLOW-CTD')
ctd_files = args_from_file(wf_dir + '/IN-FILESTOSTAGE')

data_path = '%s/data/' % wf_dir
result_path = '%s/result/' % wf_dir
log_path = '%s/logs/' % wf_dir
db_path = os.path.join(wf_dir, 'ref')

mzmlFiles = []

# mzml files
for filePath in ctd_files['Mass Spectrometry Data']:
    fileName = filePath.split('/')[-1]
    mzmlFiles.append('%s%s' % (data_path, fileName))

# Parameters
fmt = float(ctd_params['fmt'])
pmt = float(ctd_params['pmt'])
fbo = float(ctd_params['fbo'])
fdr = float(ctd_params['fdr'])
maxmod = int(ctd_params['maxmod'])

num_hits = int(ctd_params['noh'])
dmr = ctd_params['dmr']
msLevels = ctd_params['ms_levels']
centroided = ctd_params['centroided']

fixed_mods = []
variable_mods = []

mod_dict = {'fixed_mod1': 'Carbamidomethyl (C)', 'variable_mod1': 'Oxidation (M)', 'variable_mod2': 'Phospho (S)', 'variable_mod3': 'Phospho (T)', 'variable_mod4': 'Phospho (Y)'}

for c in ctd_params:
    if 'fixed' in c:
        if ctd_params[c] == 'true':
            fixed_mods.append(mod_dict[c])
    elif 'variable' in c:
        if ctd_params[c] == 'true':
            variable_mods.append(mod_dict[c])

fixed = ' '.join(fixed_mods)
variable = ' '.join(variable_mods)

logfilename = 'ligandomicsID_v2_1_coPro_workflow.logs'
logfile = open(logfilename, 'w')

if ctd_files['db'] != '':
    fasta_path = os.path.join(db_path, ctd_files['db'].split('/')[-1])
else:
    fasta_path = os.path.join(data_path, ctd_files['Individualized Reference'].split('/')[-1])

fasta_decoy_path = os.path.join(data_path, 'proteinswithdecoys.fasta')

############################################



#generate reversed decoy database
commandDecoy = 'DecoyDatabase  -in {i} -out {o} -decoy_string DECOY -decoy_string_position prefix'.format(i=fasta_path, o=fasta_decoy_path)
subprocess.call(commandDecoy.split(), stderr=logfile, stdout=logfile)

idFiles = []

#for loop over replicate mzmls
for mzml in mzmlFiles:
    if mzml.endswith('.gz'):
        logfile.write("Extracting gzipped content... \n")
        cmd = "gzip -d {f}".format(f=mzml)
        os.system(cmd)
        mzml = mzml.replace('.gz', '')
    idPath = mzml.replace('mzML', 'idXML')
    identifier = mzml.split('/')[-1].split('.')[0]

    #run peak picker if not centroided
    if centroided == 'false':
        pickpeakcommand = 'PeakPickerHiRes -in {i} -out {o} -threads 5 -algorithm:ms_levels {ms}'
        subprocess.call(pickpeakcommand.format(i=mzml, o=mzml,ms=msLevels).split(),stderr=logfile, stdout=logfile)

    #peptide search using comet
    commandComet = 'CometAdapter -in {i} -out {o} -threads 5 -database {d} -precursor_mass_tolerance {pmt} -fragment_bin_tolerance {fmt} -fragment_bin_offset {fbo} -num_hits {n} -digest_mass_range {dmr} -max_variable_mods_in_peptide {maxmod} -allowed_missed_cleavages 0'.format(i=mzml, o=idPath, d=fasta_decoy_path, pmt=pmt, fmt=fmt, fbo=fbo, n=num_hits, dmr=dmr, maxmod=maxmod) 
    if len(fixed) > 0 and len(variable) > 0:
        subprocess.call(commandComet.split() + ["-fixed_modifications", fixed, "-variable_modifications", variable, "-enzyme", "unspecific cleavage"],stderr=logfile, stdout=logfile)
    elif len(fixed) > 0:
        subprocess.call(commandComet.split() + ["-fixed_modifications", fixed, "-enzyme", "unspecific cleavage"],stderr=logfile, stdout=logfile)
    elif len(variable) > 0:
        subprocess.call(commandComet.split() + ["-variable_modifications", variable, "-enzyme", "unspecific cleavage"],stderr=logfile, stdout=logfile)
    else:
        subprocess.call(commandComet.split() + ["-enzyme", "unspecific cleavage"],stderr=logfile, stdout=logfile)

    #index decoy and target hits
    peptideIndexer = 'PeptideIndexer -in {f} -out {o} -threads 5 -fasta {d} -decoy_string DECOY -enzyme:specificity none -enzyme:name '.format(f=idPath, o=idPath, d=fasta_decoy_path)
    subprocess.call(peptideIndexer.split()  + ["unspecific cleavage"],stderr=logfile, stdout=logfile)

    #calculate fdr
    falseDiscovery = 'FalseDiscoveryRate -in {f} -out {o} -threads 5'.format(f=idPath,o=idPath.replace('.idXML','_aligner.idXML'))
    subprocess.call(falseDiscovery.split(),stderr=logfile, stdout=logfile)

    #filter by provided FDR value generate aligner maps for map transformation
    idFilter = 'IDFilter -best:strict  -in {f} -out {o} -score:pep {m} -remove_decoys -threads 5'.format(f=idPath.replace('.idXML','_aligner.idXML'), o=idPath.replace('.idXML','_aligner.idXML'), m=fdr)
    subprocess.call(idFilter.split(),stderr=logfile, stdout=logfile)
    idFiles.append(idPath)

    #id based map transformation genereate aligned idXML maps and trafoXMLs
    alignerMaps=[a.replace('.idXML', '_aligner.idXML') for a in idFiles]
    alignedMaps=[a.replace('.idXML', '_aligned.idXML') for a in idFiles]
    trafoMaps=[a.replace('.idXML', '.trafoXML') for a in idFiles]
    MapAlignerID = 'MapAlignerIdentification -in {f} -out {o} -trafo_out {t} -threads 5'.format(f=" ".join(alignerMaps), o=" ".join(alignedMaps), t=" ".join(trafoMaps))
    subprocess.call(MapAlignerID.split(), stderr=logfile, stdout=logfile)


#align mzml using trafoXMLs
for Map in trafoMaps:
    mzml = mzmlFiles[trafoMaps.index(Map)]
    trafoMZML=mzml #.replace('.mzML', '_trafo.mzML')
    MapTransformer = 'MapRTTransformer -in {f} -out {o} -trafo_in {m} -threads 5'.format(m=Map, f=mzml, o=trafoMZML)
    subprocess.call(MapTransformer.split(), stderr=logfile, stdout=logfile)

#align all content idXMLs using trafoXMLs
for Map in trafoMaps:
    ID = mzmlFiles[trafoMaps.index(Map)].replace('.mzML', '.idXML')
    trafoID = ID
    MapTransformer = 'MapRTTransformer -in {f} -out {o} -trafo_in {m} -threads 5'.format(m=Map, f=ID, o=trafoID)
    subprocess.call(MapTransformer.split(), stderr=logfile, stdout=logfile)


#merge all ids
identifer_for_files = idFiles[0].split('/')[-1]    
idresult_all = os.path.join(result_path, identifer_for_files.replace('.idXML','_merged_all.idXML'))
FileMerger = 'IDMerger -in {f} -out {o} -threads 5 -annotate_file_origin'.format(f=' '.join(idFiles),o=idresult_all)
subprocess.call(FileMerger.split(),stderr=logfile, stdout=logfile)

if num_hits >= 2:
    #filter for rank2 hits
    idresult_rank2 = os.path.join(result_path, identifer_for_files.replace('.idXML','_merged_rank2.idXML'))
    idFilter = 'IDFilter -best:n_to_m_peptide_hits 2:3  -in {f} -out {o} -score:pep 999 -threads 5'.format(f=idresult_all, o=idresult_rank2)
    subprocess.call(idFilter.split(),stderr=logfile, stdout=logfile)

    #export merged rank2 ids and all metavalues
    convert = 'TextExporter -in {f} -out {o} -id:add_hit_metavalues 0 -id:add_metavalues 0 -id:peptides_only'.format(f=idresult_rank2, o=idresult_rank2.replace('.idXML', '.csv'))
    subprocess.call(convert.split(),stderr=logfile, stdout=logfile)

#filter for rank1 hits
idresult = os.path.join(result_path, identifer_for_files.replace('.idXML','_merged.idXML'))
idFilter = 'IDFilter -best:strict  -in {f} -out {o} -score:pep 999 -threads 5'.format(f=idresult_all, o=idresult)
subprocess.call(idFilter.split(),stderr=logfile, stdout=logfile)

#export merged rank1 ids and all metavalues
convert = 'TextExporter -in {f} -out {o} -id:add_hit_metavalues 0 -id:add_metavalues 0 -id:peptides_only'.format(f=idresult, o=idresult.replace('.idXML', '_rank1.csv'))
subprocess.call(convert.split(),stderr=logfile, stdout=logfile)   

#calculate fdr on merged rank1 hits
idresult_fdr = os.path.join(result_path, identifer_for_files.replace('.idXML', '_merged_fdr.idXML'))
falseDiscovery = 'FalseDiscoveryRate -in {f} -out {o} -threads 5 -algorithm:add_decoy_peptides'.format(f=idresult,o=idresult_fdr)
subprocess.call(falseDiscovery.split(),stderr=logfile, stdout=logfile)

#extract Percolator Features with PSMFeatureExtractor
idresult_fdr_psm = os.path.join(result_path, identifer_for_files.replace('.idXML', '_merged_fdr_psm.idXML'))
PSMFeat = 'PSMFeatureExtractor -in {f} -out {o} -threads 5'.format(f=idresult_fdr,o=idresult_fdr_psm)
subprocess.call(PSMFeat.split(),stderr=logfile, stdout=logfile)

#run Percolator
idresult_perc = os.path.join(result_path, '{}_merged_perc.idXML'.format(identifer_for_files))
Percolator = 'PercolatorAdapter -in {f} -out {o} -decoy-pattern DECOY -debug 10 -threads 5 -enzyme no_enzyme -trainFDR 0.05 -testFDR 0.05'.format(f=idresult_fdr_psm,o=idresult_perc)
subprocess.call(Percolator.split(),stderr=logfile, stdout=logfile)

#filter by provided FDR value
idresult_filtered = os.path.join(result_path, '{}_merged_perc_fdr_filtered.idXML'.format(identifer_for_files))
idFilter = 'IDFilter  -in {f} -out {o} -score:pep {m} -remove_decoys -threads 5'.format(f=idresult_perc, o=idresult_filtered, m=fdr)
subprocess.call(idFilter.split(),stderr=logfile, stdout=logfile)

#export merged filtered ids and all metavalues
convert = 'TextExporter -in {f} -out {o} -id:add_hit_metavalues 0 -id:add_metavalues 0 -id:peptides_only'.format(f=idresult_filtered, o=idresult_filtered.replace('.idXML', '.csv'))
subprocess.call(convert.split(),stderr=logfile, stdout=logfile)

#find features for merged filtered ids using feature finder identification
features = []
files = [i for i in idFiles]

#assume all ids as internal ids
for file in files:
    internal=idresult_filtered
    mergeresult = os.path.join(result_path, file.replace('idXML', 'featureXML').split('/')[-1])
    features.append(mergeresult)
    IDMapper = 'FeatureFinderIdentification -in {f} -id {i} -threads 5 -out {o}'.format(f=file.replace('.idXML','.mzML'),i=internal,o=mergeresult)
    subprocess.call(IDMapper.split(),stderr=logfile, stdout=logfile)

#Link features of all replicate mzmls to one consensus file
mergeresult = os.path.join(result_path, idFiles[0].split('/')[-1].replace('idXML', 'consensusXML'))
FeatLinker = 'FeatureLinkerUnlabeledKD -in {f} -out {o} -threads 5'.format(f=' '.join(features), o=mergeresult)
subprocess.call(FeatLinker.split(),stderr=logfile, stdout=logfile)

#Resolve conflicts of different ids matching the same feature
ConfSolver = 'IDConflictResolver -in {f} -out {o}'.format(f=mergeresult,o=mergeresult)
subprocess.call(ConfSolver.split(),stderr=logfile, stdout=logfile)

#export consensus feature content
convert = 'TextExporter -in {f} -out {o} -id:add_hit_metavalues 0'.format(f=mergeresult, o=mergeresult.replace('.consensusXML', '.csv'))
subprocess.call(convert.split(),stderr=logfile, stdout=logfile)

#format output csvs
op=open(mergeresult.replace('.consensusXML', '.csv'))
opr=op.readlines()
op.close()

df = []

for i, r in enumerate(opr):
    if r.startswith('#PEPTIDE'):
        header = r.strip().split('\t')[1:] + opr[i - 1].strip().split('\t')[1:]
    if r.startswith('PEPTIDE'):
        if not opr[i-1].startswith('PEPTIDE'):
            df.append(r.strip().split('\t')[1:] + opr[i - 1].strip().split('\t')[1:])

df=pd.DataFrame(df)
df.columns=header
df.to_csv(mergeresult.replace('.consensusXML', '_edit.csv'))

#join consensus feature ids with merged rank1 ids on sequence, get number of psms and write all in one table
idexport_filtered=pd.read_csv(idresult_filtered.replace('.idXML', '.csv'), sep='\t')
idexport_rank1=pd.read_csv(idresult.replace('.idXML', '_rank1.csv'), sep='\t')
idexport_edit=pd.read_csv(mergeresult.replace('.consensusXML', '_edit.csv'), sep=',')
idexport_edit=idexport_edit[['sequence', 'rt_cf', 'mz_cf', 'intensity_cf']]
idexport_filtered=idexport_filtered.loc[idexport_filtered.groupby('sequence')['score'].idxmin()]
idexport_filtered=idexport_filtered[['sequence','score','COMET:deltCn','expect_score','file_origin','spectrum_reference','COMET:IonFrac','MS:1002252']]
merged=idexport_edit.merge(idexport_filtered, on='sequence', how='left')
num=idexport_rank1.groupby('sequence').apply(len).reset_index()
num.columns=['sequence','num_psms']
merged=merged.merge(num, on='sequence', how='left')
idexport_rank1=idexport_rank1.loc[idexport_rank1.groupby('sequence')['score'].idxmin()]
idexport_rank1=idexport_rank1[['sequence','accessions']]
merged=merged.merge(idexport_rank1, on='sequence', how='left')
if num_hits>=2:
    idexport_rank2=pd.read_csv(idresult.replace('.idXML', '_rank2.csv'), sep='\t')
    idexport_rank2=idexport_rank2[['sequence','accessions','score','MS:1002252','spectrum_reference','file_origin']]
    idexport_rank2.columns=['sequence_rank2','accessions_rank2','score_rank2','MS:1002252_rank2','spectrum_reference','file_origin']
    merged=merged.merge(idexport_rank2, on=['spectrum_reference','file_origin'], how='left')
    merged=merged[['spectrum_reference','file_origin','sequence', 'score', 'rt_cf', 'mz_cf', 'intensity_cf','accessions','num_psms','expect_score','COMET:IonFrac','MS:1002252','sequence_rank2','accessions_rank2','MS:1002252_rank2']]
    merged.columns=['spectrum_reference','file_origin','sequence', 'fdr', 'rt_cf', 'mz_cf', 'intensity','accessions','num_psms','expect_score','COMET:IonFrac','XCorr','sequence_rank2','accessions_rank2','XCorr_rank2']
    merged['deltaCn']=(merged['XCorr']-merged['XCorr_rank2'])/merged['XCorr']
else:
    merged=merged[['spectrum_reference','file_origin','sequence', 'score', 'rt_cf', 'mz_cf', 'intensity_cf','accessions','num_psms','expect_score','COMET:IonFrac','MS:1002252']]
    merged.columns=['spectrum_reference','file_origin','sequence', 'fdr', 'rt_cf', 'mz_cf', 'intensity','accessions','num_psms','expect_score','COMET:IonFrac','XCorr']
merged.to_csv(mergeresult.replace('.consensusXML', '_final_output.csv'))

logfile.close()
