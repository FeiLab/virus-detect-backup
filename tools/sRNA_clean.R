#trimRead实际中不用,测试看结果专用
trimRead <- function(fastqfile, trimmed_file, null_file, trimmed3End_file, unmatched_file, mismatch= 0.1, PCR2rc)
{

	reads <- readFastq(fastqfile);#读入FASTQ文件
	seqs <- sread(reads); # 取出所有记录（read）的sequence信息,数据格式DNAStringSet	
	rawLength <- max(width(reads));# 得到原始数据中reads长度
	lineofresult <- c(fastqfile, rawLength, length(reads))#收集3列信息,样本文件名称、reads长度和reads总数
	
	#去掉reads中含有的部分或整体PCR2rc/adapter
	trimmedCoords <- trimLRPatterns(Rpattern = PCR2rc, subject = seqs, max.Rmismatch= mismatch, with.Rindels=T,ranges=T)#这里得到坐标
	trimmedReads <- narrow(reads, start=start(trimmedCoords), end=end(trimmedCoords))#利用上一步得到的坐标，同时trim核苷酸序列和质量分数序列	
	trimmed3End <- narrow(reads, start=end(trimmedCoords)+1, end=width(reads))#把trimm掉的那部分序列保留，以备人工检查
	rm(trimmedCoords)
	writeFastq(trimmed3End, file=trimmed3End_file)      	
	rm(trimmed3End)
	
	unmatchedReads <- trimmedReads[width(trimmedReads)==width(reads)[1]]#保存全长的reads，表示unmatched
	writeFastq(unmatchedReads, file= unmatched_file)
	unmatched_length <- length(unmatchedReads)
	rm(unmatchedReads)
	nullReads <- trimmedReads[width(trimmedReads)==0]#保存空reads，表示null
	writeFastq(nullReads, file= null_file)
	null_length <- length(nullReads)
	rm(nullReads)
	trimmedReads <- trimmedReads[width(trimmedReads)<width(reads)[1]]#保存非全长的reads，表示trimmed
	trimmedReads <- trimmedReads[width(trimmedReads)>0]#去掉无序列的空数据
	writeFastq(trimmedReads, file=trimmed_file)
	lineofresult <- c(lineofresult, length(trimmedReads),null_length,unmatched_length)#加入1列（trimmedReads总数）
	write(lineofresult,file = "trimmed.report", ncolumns =7,append = T, sep = "\t")
	rm(trimmedReads)
	gc()
}

cleanRead <- function(fastqfile, cleaned_file, nCutoff=1, readLength=15, RdPerYield){
	#reads <- readFastq(fastqfile);#读入FASTQ文件
	inFh <- FastqStreamer(fastqfile, n=RdPerYield); 
	if (file.exists(cleaned_file) ) {file.remove(cleaned_file); } #如果输出文件已经存在必须删除，防止追加写
	iteration=0;
	trimmed_reads=0;
	cleaned_reads=0;
	trimmed_len=0;
	cleaned_len=0;
	while (length(reads <- yield(inFh))) { #每次控制读入5百万个reads
		seqs <- sread(reads); # 取出所有记录（read）的sequence信息,数据格式DNAStringSet		
		nCount<-alphabetFrequency(seqs)[,"N"];#统计每条read中的字符N总数,
		rm(seqs);#用完及时删除		
		trimmed_reads <- trimmed_reads + length(reads);#累计trimmed的reads的总数
		trimmed_len=trimmed_len+sum(width(reads));#累计trimmed reads的总长度		
		cleanedReads<-reads[nCount<nCutoff];#只保留字符N总数<nCutoff的reads
		rm(reads);
		cleanedReads <- cleanedReads[width(cleanedReads)>=readLength];#去掉不足一定长度（默认是15bp）的reads
		writeFastq(cleanedReads, file=cleaned_file, mode="a", full=FALSE);#保存cleaned的数据,注意是追加写
		cleaned_reads <- cleaned_reads + length(cleanedReads);#累计cleaned的reads的总数
		cleaned_len=cleaned_len+sum(width(cleanedReads));#累计cleaned reads的总长度
		rm(cleanedReads);
		gc();
	}#End yield while;
	close(inFh);
	trimmed_len=trimmed_len/trimmed_reads;#得到trimmed reads的平均长度
	cleaned_len=cleaned_len/cleaned_reads;#得到cleaned reads的平均长度
	lineofresult <- c(fastqfile, trimmed_reads, trimmed_len, cleaned_reads, cleaned_len);
	write(lineofresult,file = "cleaned.report", ncolumns =5,append = T, sep = "\t");
}
