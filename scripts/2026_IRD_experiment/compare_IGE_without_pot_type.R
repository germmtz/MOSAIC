# Sensitivity of M3 predictions to omission of the mono/mixt FIXED effect.
# Run: Rscript scripts/2026_IRD_experiment/compare_IGE_without_pot_type.R
# Refit M3 for traits selected as M3 in the primary analysis; do not reselect models.
# All primary results stay unchanged. Solo predictions are not refitted: their
# model already has no pot_type contrast because all its observations are solo.
invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
root <- dirname(dirname(dirname(normalizePath(sub('^--file=','',script_arg[1])))))
primary <- file.path(root,'outputs','2026_IRD_expermiment','ige_estimations')
out <- file.path(primary,'without_pot_type')
dir.create(out,showWarnings=FALSE)
primary_files <- list.files(primary,full.names=TRUE)
primary_files <- primary_files[!dir.exists(primary_files)]
original_hashes <- tools::md5sum(primary_files)
read_csv <- function(path) read.csv2(path,fileEncoding='UTF-8-BOM',na.strings='NA')
write_csv <- function(d,name){
 lines<-character();tc<-textConnection('lines','w',local=TRUE)
 write.csv2(d,tc,row.names=FALSE,na='NA');close(tc)
 con<-file(file.path(out,name),'wb');on.exit(close(con))
 writeBin(as.raw(c(239,187,191)),con)
 writeBin(charToRaw(enc2utf8(paste0(paste(lines,collapse='\r\n'),'\r\n'))),con)
}
selection <- read_csv(file.path(primary,'01_selected_models.csv'))
traits <- selection$trait[selection$selected_model=='M3']
if(!length(traits))stop('No primary M3 traits to compare.')
primary_dge <- read_csv(file.path(primary,'03_DGE_genotype_predictions.csv'))
primary_pairs <- read_csv(file.path(primary,'03_all_genotype_interactions.csv'))
all_dge<-all_pairs<-all_vc<-all_starts<-all_diagnostics<-list()
fit_once <- function(d,gp=NULL,rp=NULL){
 asreml::asreml(fixed=response~chapelle+resown,
  random=~row_group+column_group+str(~DGE+IGE,~corgh(2):id(DGE))+DGExIGE,
  residual=~id(pot):corv(plant),data=d,
  G.param=gp,R.param=rp,
  na.action=asreml::na.method(y='include',x='fail'),
  maxit=100,workspace='256mb',trace=FALSE)
}
for(t in traits){
 message('Fitting M3 without pot_type: ',t)
 d<-read_csv(file.path(primary,paste0('input_',t,'.csv')))
 d<-d[order(d$pot_id,d$plant_id),]
 genotypes<-sort(unique(c(d$genotype_id,d$neighbour_id)))
 d$DGE<-factor(d$genotype_id,levels=genotypes)
 d$IGE<-factor(d$neighbour_id,levels=genotypes)
 d$DGExIGE<-factor(paste0('g',d$genotype_id,'_n',d$neighbour_id))
 d$pot<-factor(d$pot_id);d$plant<-factor(d$plant_id,levels=1:2)
 d$chapelle<-factor(d$chapelle);d$resown<-factor(d$resown)
 d$row_group<-interaction(d$chapelle,d$row,drop=TRUE)
 d$column_group<-interaction(d$chapelle,d$column,drop=TRUE)
 stopifnot(!anyDuplicated(d[c('pot','plant')]),all(table(d$pot)==2))
 fit<-fit_once(d)
 if(!isTRUE(fit$converge))fit<-update(fit,maxit=100)
 if(!isTRUE(fit$converge))stop('Default fit failed for ',t)
 starts<-data.frame(trait=t,start='default',converged=TRUE,loglik=fit$loglik)
 base<-fit
 # Match the primary analysis: negative, zero, positive correlation starts.
 for(rho in c(-.8,0,.8)){
  gp<-base$G.param
  for(term in names(gp))for(j in seq_along(gp[[term]])){
   b<-gp[[term]][[j]];if(!is.null(b$con))b$con[b$con=='B']<-'P';gp[[term]][[j]]<-b
  }
  gp[['DGE+IGE']][['2']]$initial[]<-c(rho,.15*var(d$response,na.rm=TRUE),.05*var(d$response,na.rm=TRUE))
  gp[['DGE+IGE']][['2']]$con<-c('U','P','P')
  candidate<-try(fit_once(d,gp,base$R.param),silent=TRUE)
  valid<-!inherits(candidate,'try-error')&&isTRUE(candidate$converge)
  starts<-rbind(starts,data.frame(trait=t,start=paste0('rho=',rho),converged=valid,
   loglik=if(valid)candidate$loglik else NA_real_))
  if(valid&&candidate$loglik>fit$loglik)fit<-candidate
 }
 vc<-summary(fit)$varcomp
 if('bound'%in%names(vc)&&any(as.character(vc$bound)%in%c('S','C')))stop('Singular variance block: ',t)
 saveRDS(fit,file.path(out,paste0(t,'_M3_without_pot_type.rds')))
 all_starts[[t]]<-starts
 all_vc[[t]]<-data.frame(trait=t,term=rownames(vc),vc,row.names=NULL)
 fc<-summary(fit,coef=TRUE)$coef.fixed
 write_csv(data.frame(term=rownames(fc),fc,row.names=NULL),paste0(t,'_fixed_effects.csv'))
 all_diagnostics[[t]]<-data.frame(trait=t,converged=fit$converge,loglik=fit$loglik,
   observed=sum(!is.na(d$response)),pots=nlevels(d$pot))
 getb<-function(term,level){
  b<-fit$coefficients$random;ix<-match(paste0(term,'_',level),rownames(b))
  if(anyNA(ix))stop('Unmatched BLUP: ',term);as.numeric(b[ix])
 }
 ix<-match('(Intercept)',rownames(fc));if(is.na(ix))stop('Missing intercept')
 intercept<-fc[ix,'solution'];blup<-getb('DGE',as.character(genotypes))
 all_dge[[t]]<-data.frame(trait=t,genotype=as.character(genotypes),intercept=intercept,
  DGE_BLUP=blup,prediction=intercept+blup,reference_chapelle=levels(d$chapelle)[1],
  reference_resown=levels(d$resown)[1],pot_type_fixed_effect='omitted')
 pairs<-expand.grid(genotype=as.character(genotypes),neighbour=as.character(genotypes),stringsAsFactors=FALSE)
 effect<-getb('IGE',pairs$neighbour)+getb('DGExIGE',paste0('g',pairs$genotype,'_n',pairs$neighbour))
 all_pairs[[t]]<-data.frame(trait=t,pairs,interaction=effect,interaction_type='IGE + DGE:IGE')
}
dge<-do.call(rbind,all_dge);pairs<-do.call(rbind,all_pairs)
write_csv(dge,'03_DGE_genotype_predictions.csv')
write_csv(pairs,'03_all_genotype_interactions.csv')
write_csv(do.call(rbind,all_vc),'variance_parameters.csv')
write_csv(do.call(rbind,all_starts),'starting_value_checks.csv')
write_csv(do.call(rbind,all_diagnostics),'fit_diagnostics.csv')
# Exact genotype/pair joins; do not match by row order.
stopifnot(!anyDuplicated(primary_dge[c('trait','genotype')]),
 !anyDuplicated(primary_pairs[c('trait','genotype','neighbour')]))
a<-merge(primary_dge[c('trait','genotype','prediction')],dge[c('trait','genotype','prediction')],
 by=c('trait','genotype'),suffixes=c('_with_pot_type','_without_pot_type'))
b<-merge(primary_pairs[c('trait','genotype','neighbour','interaction')],pairs[c('trait','genotype','neighbour','interaction')],
 by=c('trait','genotype','neighbour'),suffixes=c('_with_pot_type','_without_pot_type'))
stopifnot(nrow(a)==10L*length(traits),nrow(b)==100L*length(traits))
write_csv(a,'DGE_matched_comparison.csv');write_csv(b,'interactions_matched_comparison.csv')
comparisons<-list(
 DGE=data.frame(trait=a$trait,genotype=a$genotype,neighbour=NA_character_,
  original=a$prediction_with_pot_type,without=a$prediction_without_pot_type),
 same_genotype=data.frame(trait=b$trait[b$genotype==b$neighbour],genotype=b$genotype[b$genotype==b$neighbour],
  neighbour=b$neighbour[b$genotype==b$neighbour],original=b$interaction_with_pot_type[b$genotype==b$neighbour],
  without=b$interaction_without_pot_type[b$genotype==b$neighbour]),
 different_genotype=data.frame(trait=b$trait[b$genotype!=b$neighbour],genotype=b$genotype[b$genotype!=b$neighbour],
  neighbour=b$neighbour[b$genotype!=b$neighbour],original=b$interaction_with_pot_type[b$genotype!=b$neighbour],
  without=b$interaction_without_pot_type[b$genotype!=b$neighbour]))
results<-list()
for(type in names(comparisons))for(t in traits){
 xy<-comparisons[[type]];xy<-xy[xy$trait==t & complete.cases(xy[c('original','without')]),]
 test<-cor.test(xy$original,xy$without,method='pearson')
 results[[paste(type,t)]]<-data.frame(trait=t,prediction_type=type,n=nrow(xy),pearson_R=unname(test$estimate),
  p_nominal=test$p.value,mean_change=mean(xy$without-xy$original),RMSE_change=sqrt(mean((xy$without-xy$original)^2)),
  note='Descriptive; estimated effects have shared uncertainty, particularly directed pairs')
 p<-ggplot(xy,aes(original,without))+geom_abline(slope=1,intercept=0,linetype=2,colour='grey50')+
  geom_smooth(method='lm',formula=y~x,se=FALSE,colour='#B95442',linewidth=.7)+
  geom_point(colour='#277DA8',size=if(type=='different_genotype')1.8 else 2.7,alpha=.8)+
  theme_bw(base_size=12)+labs(title=paste('Omitting pot type:',gsub('_',' ',type)),
   subtitle=sprintf('Pearson R = %.3f; nominal p = %s; n = %d',test$estimate,format.pval(test$p.value,digits=3,eps=.001),nrow(xy)),
   x='Prediction with mono/mixt fixed effect (g)',y='Prediction without mono/mixt fixed effect (g)',
   caption='Dashed line: equality. Refitted M3; all other terms retained. Nominal inference ignores prediction uncertainty.')
 if(type!='different_genotype')p<-p+geom_text(aes(label=genotype),vjust=-.7,size=3.5)
 ggsave(file.path(out,paste0('comparison_',type,'_',t,'.png')),p,width=8,height=7,dpi=400,bg='white')
}
results<-do.call(rbind,results);write_csv(results,'prediction_correlations.csv')
writeLines(c('M3 refitted without pot_type for traits selected as M3 in the PRIMARY analysis; no new selection.',
 'Fixed: chapelle + resown. Random: row/column within chapelle, correlated DGE/IGE, directed interaction.',
 'Residual: id(pot):corv(plant). All observations retained; zeros and resown plants included.',
 'DGE prediction = model intercept + focal DGE BLUP; reference chapelle 1, non-resown, other random effects zero.',
 'Interaction = neighbour IGE + directed interaction; no mono contrast in the refitted model.',
 'Primary same-genotype predictions include the original mono fixed contrast; original between-genotype predictions do not.',
 'Correlations are descriptive; nominal p-values ignore dependence and estimation uncertainty.',
 'REML likelihoods of models with DIFFERENT fixed effects must not be compared directly as an LRT.',
 'Solo model already omits pot_type; it is unchanged and was not refitted.'),file.path(out,'analysis_notes.txt'))
stopifnot(identical(unname(original_hashes),unname(tools::md5sum(primary_files))))
message('Primary output files verified unchanged.')
print(results,row.names=FALSE)