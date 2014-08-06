#######
#	updates data in targetPredictor.db
#
#######

##
# TODO:
# change paths to current path, no need to alter library
##



require(stringr)		#string operations, essential
require(sqldf) 			#sqlite dataframe operations, i don't think i use it...
require(RSQLite)		#local sqlite db, essential

require(biomaRt)		#connection to ensembl through biomart 

require(targetPredictor.db)

options(sqldf.driver = "SQLite") # as per FAQ #7 force SQLite
options(gsubfn.engine = "R")



#to put in the Updater package?
#' @title Auxiliary internal TargetPredictor.db functions
#'
#' This function is not meant to be called directly by the user
#'
#' @export
#' @author Maciej Pajak \email{m.pajak@@sms.ed.ac.uk}
get_connection2 = function() {
	tp_db=paste(system.file(package="targetPredictor.db"),"/extdata/targetpredictor_db2.sqlite",sep="")
	return(RSQLite::dbConnect(dbDriver("SQLite"), dbname = tp_db))
			#on.exit(dbDisconnect(conL))
			
}


#	pathToDb = 'C:/Users/Maciej/workspace-R/targetPredictor.db/inst/extdata/targetpredictor_db.sqlite'
metadata = function(pathToDb) {
	
	tbl = data.frame(name = c('package','Db type'), value=c('targetPredictor.db','mirnaDb'))
	conL = RSQLite::dbConnect(dbDriver("SQLite"), dbname = pathToDb)
	if (dbExistsTable(conL, 'metadata')) {
		dbRemoveTable(conL, 'metadata')
	}
	dbWriteTable(conL, 'metadata', tbl, row.names=FALSE)
	dbDisconnect(conL)
	

}



#
#
#UPDATE
#TARGETSCAN
#
#to be run on conL2 = get_connection2()
update_targetscan = function(conL,species='mmu', ...) {
	
	if (species=='mmu') {
		tax.num='10090'
		ts_url = 'http://www.targetscan.org/mmu_61/mmu_61_data_download/Conserved_Site_Context_Scores.txt.zip'
	} else if (species=='hsa') {
		tax.num='9606'
		ts_url = 'http://www.targetscan.org/vert_61/vert_61_data_download/Conserved_Site_Context_Scores.txt.zip'
	}
	
	#testing workspace
	#targetscan_dest_context = 'Conserved_Site_Context_Scores.txt'
	targetscan_dest_context.zip=paste(path.package("targetPredictor.db"),"/extdata/targetscan_context_",species,".zip",sep="")
	targetscan_dest.path=paste(path.package("targetPredictor.db"),"/extdata",sep="")
	targetscan_dest_context=paste(path.package("targetPredictor.db"),"/extdata/Conserved_Site_Context_Scores.txt",sep="")
	
	
	if (!file.exists(targetscan_dest_context.zip)) {
		download.file(url=ts_url, destfile=targetscan_dest_context.zip)
		unzip(targetscan_dest_context.zip, overwrite=TRUE, exdir=targetscan_dest.path)
		
		cols_context = c('Gene_ID', 'Gene_Symbol', 'Transcript_ID', 'Gene_Tax_ID', 'miRNA', 'Site Type', 'UTR_start', 'UTR_end', 'three_prime_pairing', 'local_AU', 'position', 'TA', 'SPS', 'context_score', 'context_score_percentile')
		#conL2
		dbWriteTable(conL, name="targetscan_context", value=targetscan_dest_context, 
					row.names=FALSE, header=FALSE, col.names=cols_context, skip=1, sep = "\t", overwrite=TRUE)
		
		mapped_table = paste('targetscan_mapped_',species,sep='')
		
		mirna_query = paste('CREATE TABLE ',mapped_table,' AS 
		SELECT DISTINCT miRNA AS mirna, -context_score AS score, Gene_ID AS GeneID FROM targetscan_context WHERE Gene_Tax_ID=', tax.num, ' AND context_score<0', 
					sep= '')
						
		#conL2			
		xo = dbGetQuery(conL, mirna_query)	
	}	
	
}

#
#
#UPDATE
#MIRANDA
#
update_miranda = function(conL,species='mmu', ...) {
	if (species=='mmu') {
		tax.num='10090'
		miranda_url = 'ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.mus_musculus.zip'
		miranda_dest=paste(path.package("targetPredictor.db"),"/extdata/v5.txt.mus_musculus",sep="")
		ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
	} else if (species=='hsa') {
		tax.num='9606'
		miranda_url = 'ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.homo_sapiens.zip'
		miranda_dest=paste(path.package("targetPredictor.db"),"/extdata/v5.txt.homo_sapiens",sep="")
		ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
	} else if (species=='rno') {
		tax.num='10116'
		miranda_url = 'ftp://ftp.ebi.ac.uk/pub/databases/microcosm/v5/arch.v5.txt.rattus_norvegicus.zip'
		miranda_dest=paste(path.package("targetPredictor.db"),"/extdata/v5.txt.rattus_norvegicus",sep="")
		ensembl = useMart("ensembl",dataset="rnorvegicus_gene_ensembl")
	}
	
	miranda_dest.zip=paste(path.package("targetPredictor.db"),"/extdata/miranda",species,".zip",sep="")
	miranda_dest.path=paste(path.package("targetPredictor.db"),"/extdata",sep="")
		
	unmapped_table=paste('miranda_',species,sep='')
	mapped_table=paste('miranda_mapped_',species,sep='')
	
	if (!file.exists(miranda_dest)) {
		message('first run of Miranda - downloading data...')
		download.file(url=miranda_url, destfile=miranda_dest.zip)	
		unzip(miranda_dest.zip,overwrite=TRUE, exdir=miranda_dest.path)
		
		cols = c('GROUP','SEQ','METHOD','FEATURE','CHR','START','END','STRAND','PHASE','SCORE','PVALUE_OG','TRANSCRIPT_ID','EXTERNAL_NAME')
		message('first run of Miranda - converting to sqlite...')
		#first method is slow but works with any EOL symbols
		miranda.temp = read.table(file=miranda_dest, header=FALSE, col.names=cols, skip=5, sep = "\t", stringsAsFactors=FALSE)
		dbWriteTable(conL, name=unmapped_table, value=miranda.temp, overwrite=TRUE)
		rm(miranda.temp)
		
		#while this method is definitely preferred, its weakness is dependence of windows EOL symbols
		#dbWriteTable(conL, name=unmapped_table, value=miranda_dest, 
		#	row.names=FALSE, header=FALSE, col.names=cols, skip=5, sep = "\t", overwrite=TRUE)
	}	
	
	mirandagenes_query = paste('SELECT DISTINCT TRANSCRIPT_ID FROM ',unmapped_table,sep='')
	eo = dbGetQuery(conL, mirandagenes_query)
	
	#if (species=='mmu') {
	#	ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
	#}
	relevantBM=getBM(attributes=c('ensembl_transcript_id','entrezgene'),filters='ensembl_transcript_id', values=eo$TRANSCRIPT_ID ,mart=ensembl)
	relevantBM=relevantBM[complete.cases(relevantBM),]
	
	dbWriteTable(conL, name=paste('mappingsENS2GeneID_',species,sep=''), value=relevantBM, overwrite=TRUE)
	
	if	(dbExistsTable(conL, mapped_table)) {
		dbRemoveTable(conL, mapped_table)
	}
	join_query = paste('CREATE TABLE ',mapped_table,' AS 
				SELECT DISTINCT miranda.SEQ as mirna, miranda.SCORE as score, mapping.entrezgene AS GeneID 
				FROM ',unmapped_table,' AS miranda, mappingsENS2GeneID_',species,' AS mapping
				WHERE miranda.TRANSCRIPT_ID=mapping.ensembl_transcript_id',sep='')
				
	eo = dbGetQuery(conL, join_query)
	
	if	(dbExistsTable(conL, mapped_table) && dbExistsTable(conL, unmapped_table)) {
		dbRemoveTable(conL, unmapped_table)
	}
	
	return(0)
}


update_diana_1s = function(full_data, conL, species='mmu', ... ) {
	
	print(species)
	
	query = paste('SELECT Gene, miR, score 
			FROM full_data
			WHERE score>0.3 AND miR LIKE "',species,'%"',sep='')
			
	eo = sqldf(query)
	
	print('species extracted')
	
	eo$Gene = str_match(as.character(eo$Gene),'([A-Z0-9]*?)\\(.*?\\)')[,2]
	print('gene names extracted')
	
	eo$miR = str_match(as.character(eo$miR),paste('(',species,			'-(?:miR|let)-[a-z0-9]*?(?:-[a-z0-9]*?)?)\\(.*?\\)',sep=''))[,2]
	print('mirna names extracted')
	
	unmapped_table = paste('diana_',species,sep='')
	mapped_table = paste('diana_mapped_',species,sep='')
	
	dbWriteTable(conL, name=unmapped_table,
				value=eo, 
				row.names=FALSE, overwrite=TRUE)
	print('unmapped table written, mapping names...')
	
	if	(dbExistsTable(conL, mapped_table)) {
		dbRemoveTable(conL, mapped_table)
	}			
	join_query = paste('CREATE TABLE ',mapped_table,' AS 
				SELECT DISTINCT diana.miR AS mirna, diana.score AS score, mapping.GeneID 
				FROM ',unmapped_table,' AS diana, mappingsRNA2GeneID_',species,' AS mapping
				WHERE diana.Gene=mapping.Ensembl_gene_identifier',sep='')
				
	eo = dbGetQuery(conL, join_query)				
}


#
#
#UPDATE
#DIANA
#
update_diana_all = function(conL, ...) {
#from a flat file, i cannot believe it i haven't found it earlier
	diana_url = 'http://diana.imis.athena-innovation.gr/DianaTools/data/microtv4_data.tar.gz'
	
	WD = "C:/Users/Maciej/workspace-R"
	setwd(WD)
	fullcsv='microtv4_data.csv'
	
	xo = read.table(file=fullcsv,fill=TRUE, sep=',', header=TRUE)
	#that's an abuse but after half an hour it reads the data frame in
	xo=xo[complete.cases(xo),]
	
	colnames(xo) = c('Transcipt','Gene','miR','score')
	
	update_diana_1s(xo,conL,species='mmu')
	update_diana_1s(xo,conL,species='hsa')
	
	return(0)
}




#
#
#UPDATE
#PICTAR
#
update_pictar = function(conL, species='mmu', ...) {

		#
	if (species=='mmu') {
		tax.num='10090'
		pictar_url = 'http://dorina.mdc-berlin.de/rbp_browser/downloads/pictar_mm9_mammals.bulk_download.csv'
	} else if (species=='hsa') {
		tax.num='9606'
		pictar_url = 'http://dorina.mdc-berlin.de/rbp_browser/downloads/	pictar_hg19_mammals.bulk_download.csv'
	}
	#tax.rno='10116'

	pictar_dest=paste(path.package("targetPredictor.db"),"/extdata/pictar",species,".csv",sep="")	
	#pictar_temp = paste('C:/Users/Maciej/workspace-R/pictar',species,'.csv',sep='')
	
	unmapped_table=paste('pictar_',species,sep='')
	mapped_table=paste('pictar_mapped_',species,sep='')

	if (!file.exists(pictar_dest)) {
			message('first run of Pictar - downloading data...')
			download.file(url=pictar_url, destfile=pictar_dest)
			cols=c('target','score','mirna')
			# read csv file into sql database (good way of dealing with huge files without engaging RAM)
			
			message('first run of Pictar - converting to sqlite...')
			dbWriteTable(conL, name=unmapped_table,
				value=pictar_dest, 
				row.names=FALSE, header=FALSE, col.names=cols, skip=1, sep = ",", overwrite=TRUE)
	}
		
	
	#
	#	JOINING PICTAR:

	query = paste("SELECT GeneID, RNA_nucleotide_accession_version, Ensembl_rna_identifier, Ensembl_gene_identifier FROM gene2ensembl WHERE tax_id='", 
					tax.num,"'",sep='')
	
	relevant = dbGetQuery(conL, query)
	relevant$RNArefseq = str_match(relevant$RNA_nucleotide_accession_version, '[A-Z_0-9]+')
	dbWriteTable(conL, name=paste('mappingsRNA2GeneID_',species,sep=''),
			value=relevant[,c('GeneID','RNArefseq','Ensembl_rna_identifier','Ensembl_gene_identifier')], overwrite=TRUE)
	

	
	if	(dbExistsTable(conL, mapped_table)) {
		dbRemoveTable(conL, mapped_table)
	}
	join_query = paste('CREATE TABLE ',mapped_table,' AS 
				SELECT DISTINCT pictar.mirna AS mirna, pictar.score AS score, mapping.GeneID 
				FROM ',unmapped_table,' AS pictar, mappingsRNA2GeneID_',species,' AS mapping
				WHERE pictar.target=mapping.RNArefseq',sep='')
				
	eo = dbGetQuery(conL, join_query)

	if	(dbExistsTable(conL, mapped_table) && dbExistsTable(conL, mapped_table)) {
		dbRemoveTable(conL, unmapped_table)
	}
	
	return(0)
}



#
#
#	UPDATE
#	CLASH
#	(from Cell Hewlak et al. 2013)
#
#
update_clash = function(conL, study='hewlak', species='hsa',...) { 
	#clash doesn't need updating as such, the name is just for consistency
	#study var might be useful at some point
	species='hsa'	#it is hsa bcoz it's just one study, but what if we have access to more in other species...
	tax.num='9606'
		
	cols = c('seq_ID', 'microRNA_name', 'miRNA_start', 'miRNA_end', 'miRNA_seq', 'mRNA_name', 'mRNA_start', 'mRNA_end_extended', 'mRNA_seq_extended', 'chimeras_decompressed', 'experiments', 'experiments_list', 'microRNA_first', 'two_way_merged', 'seed_type', 'num_basepairs', 'seed_basepairs', 'folding_energy', '5_UTR', 'CDS', '3_UTR', 'folding_class', 'conservation_score', 'log2_target_enrichment', 'CLASH_single_reads_ovlp', 'CLASH_cluster_ovlp', 'PAR_CLIP_cluster_ovlp')	
		
	clash.data = read.table(file='clash-interactions.txt',sep='\t', skip=31, header=FALSE, col.names=cols, stringsAsFactors=FALSE)
	#temporary reading local file

	clash.data$miR = paste(species,'-',str_match(as.character(clash.data$microRNA_name),
								'(?:.*?)\\_((?:let|miR)-[a-z0-9]+)')[,2],sep='')
	clash.data$geneENS = str_match(as.character(clash.data$mRNA_name),
								'(ENSG[0-9]+)')[,2]
	clash.slim = clash.data[,c('miR','geneENS','folding_energy','folding_class')]
	write.csv(clash.slim, file='clash_slim.csv')

	dbWriteTable(conL, name='clash_hewlak_hsa', value=clash.slim, overwrite=TRUE)

	#join 

	query = paste("SELECT GeneID, RNA_nucleotide_accession_version, Ensembl_rna_identifier, 	
					Ensembl_gene_identifier FROM gene2ensembl WHERE tax_id='", 
						tax.num,"'",sep='')
	
	relevant = dbGetQuery(conL, query)
	relevant$RNArefseq = str_match(relevant$RNA_nucleotide_accession_version, '[A-Z_0-9]+')
	dbWriteTable(conL, name='mappingsRNA2GeneID_hsa', value=relevant[,c('GeneID','RNArefseq','Ensembl_rna_identifier','Ensembl_gene_identifier')], overwrite=TRUE)

	if	(dbExistsTable(conL, 'clash_hewlak_mapped_hsa')) {
			dbRemoveTable(conL,'clash_hewlak_mapped_hsa')
		}
		join_query = 'CREATE TABLE clash_hewlak_mapped_hsa AS 
					SELECT DISTINCT clash.miR AS mirna, -clash.folding_energy AS score, mapping.GeneID 
					FROM clash_hewlak_hsa AS clash, mappingsRNA2GeneID_hsa AS mapping
					WHERE clash.geneENS=mapping.Ensembl_gene_identifier'
					
		eo = dbGetQuery(conL, join_query)

	#if	(dbExistsTable(conL, 'clash_hewlak_mapped_hsa') && dbExistsTable(conL, 'clash_hewlak_hsa')) {
	#		dbRemoveTable(conL,'clash_hewlak_hsa')
	#	}
	
	
}


#
#
#	UPDATE
#	MAPPINGS
#
#
update_mappings = function(conL,...) {
	if (species=='mmu') {
		tax.num='10090'
	} else if (species=='hsa') {
		tax.num='9606'
	}

	gene2ensembl_url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz' #directly read it from online to a local sqlite?
	mapping_dest.gz=paste(path.package("targetPredictor.db"),"/extdata/gene2ensembl.gz",sep="")
	mapping_dest=paste(path.package("targetPredictor.db"),"/extdata/gene2ensembl",sep="")
	
		
	if (!file.exists(mapping_dest)) {
	
		download.file(url=gene2ensembl_url, destfile=mapping_dest.gz)
		
		system(paste('gzip -d -f ',path.package("targetPredictor.db"),'/extdata/gene2ensembl.gz', sep='')) 
		#not sure about direct use of shell but R's untar and unzip actually use shell too so who cares
		
		cols = c('tax_id', 'GeneID', 'Ensembl_gene_identifier', 'RNA_nucleotide_accession_version', 'Ensembl_rna_identifier', 'protein_accession_version', 'Ensembl_protein_identifier')

		dbWriteTable(conL, name="gene2ensembl", value=mapping_dest, 
				 row.names=FALSE, header=FALSE, col.names=cols, skip=1, sep = "\t", overwrite=TRUE)
	
	} 
}

#
#
#	UPDATE
#	HOMOLOGENE
#
#
update_homologene = function(conL,...) {
	
	current=FALSE
	
	if (current==FALSE) {
		homologene_dest=paste(path.package("targetPredictor.db"),"/extdata/homologene.data",sep="")
		homologene_url = 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data'
		download.file(url=homologene_url, destfile=homologene_dest)
		
		cols=c('family','tax_id','entrez','symbol','pm','pr')
		
		dbWriteTable(conL, name="homologene_raw", value=homologene_dest, 
					row.names=FALSE, header=FALSE, col.names=cols, skip=0, sep = "\t", overwrite=TRUE)
	
	return(0)
}
}

update_mirtar = function(conL,...) { 
	#BSH, works only if the right csv is the workspace
	mirtar = read.csv(file='miRTarBase_MTI.csv', stringsAsFactors=FALSE)
	full_species=c('Drosophila melanogaster','Rattus norvegicus','Mus musculus', 'Homo sapiens')
	mirtar = mirtar[mirtar$Species..miRNA. %in% full_species,c(2,5)]
	colnames(mirtar) = c('mirna','GeneID')
	
	dbWriteTable(conL,name='mirtar_all',value=mirtar, row.names=FALSE, header=TRUE)
	

}	
	
refresh_data = function() {
	
	conL2 = get_connection2()
	
	#deletes flat files so with the next run you will get updated ones
	#update_diana(all.rnas,conL,species='mmu')
	update_mappings(conL2)
	update_homologene(conL2)
	update_targetscan(conL2, species='mmu')
	update_targetscan(conL2, species='hsa')
	
	update_miranda(conL2, species='mmu')
	update_miranda(conL2, species='hsa')
	update_miranda(conL2, species='rno') #THE ONLY ONE
	
	update_pictar(conL2, species='mmu')
	update_pictar(conL2, species='hsa')
	
	update_diana(conL2, species='mmu')
	update_diana(conL2, species='hsa')
	
	
	return(0)
}


copy_table = function(from, to, tabname) {
	dbWriteTable(to, name=tabname, value=dbReadTable(from,tabname), overwrite=TRUE, row.names=FALSE)
}

#slimming the db
ss0 = function() {

	conL = get_connection()
	conLL = get_connection2()

	interesting_tabs = c("clash_hewlak_mapped_hsa",'diana_mapped_mmu','diana_mapped_hsa',
						"mappingsENS2GeneID_hsa", "mappingsENS2GeneID_mmu", "mappingsRNA2GeneID_mmu",
						"mappingsRNA2GeneID_hsa", "homologene_raw", "pictar_mapped_mmu",
						"pictar_mapped_hsa", "miranda_mapped_mmu", "miranda_mapped_hsa",
						"miranda_mapped_rno", "targetscan_mapped_hsa",   "targetscan_mapped_mmu")
						
	for  (tab in interesting_tabs) {
		copy_table(conL,conLL,tab)
	}

}
