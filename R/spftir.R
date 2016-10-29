#' @title Truncation of a Region of a Mid-Infrared Spectral Matrix
#' @description Allow to trim a region of the spectra defined between two wavenumbers.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param limInf numeric. Upper wavenumber limit of the spectral region.
#' @param limSup numeric. Lower wavenumber limit of the spectral region.
#' @return A truncated matrix within two wavenumber limits. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Truncated
#' trn <- sptrun(spectra=spectra, limInf=800, limSup=2000)
#' @export

sptrun <- function(spectra, limInf, limSup){

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(limSup)) {
    stop('No upper spectral limit provided')
  }
  if (missing(limInf)) {
    stop('No lower spectral limit provided')
  }
  if (is.null(which(round(spectra[1, ], 0) == limSup)) == TRUE) {
    stop('No upper limit spectral provided')
  }
  if (is.null(which(round(spectra[1, ], 0) == limInf)) == TRUE) {
    stop('No lower limit spectral provided')
  }

  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }

  LS <- which(round(spectra[1, ], 0) == limSup)[1]
  LI <- which(round(spectra[1, ], 0) == limInf)[1]
  output <- spectra[ , LS:LI]
}


#' @title Interpolation for Intermediate Values of a Matrix of Mid-infrared Spectra
#' @description Allow the interpolation of intermediate values of a matrix of mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param n numeric. Number of interpolated values between two variables. Defaults to 1.
#' @importFrom stats approx
#' @return A matrix spectra with interpolated values. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Interpolated spectra
#' int <- spint(spectra=spectra, n=1)
#' @export

spint<-function(spectra, n=1){

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(n)) {
    stop('No n value provided')
  }

  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  k<-(length(spectra[1,])-1)*n + length(spectra[1,])
  s<-matrix(NA, nrow=2, ncol=k)
  s[1,1]<-spectra[1,1]
  a=1
  for (i in 1:(length(spectra[1,])-1)){
    a<-a+n+1
    s[1,a]<- spectra[1,i+1]
    for(j in 1:n){
      s[1,(a-j)]<-spectra[1,i+1]-(spectra[1,i+1]-spectra[1,i])/(n+1)*j
    }
  }
  out<-approx(spectra[1,] , spectra[2,] , xout = s[1,], method="linear")
  output <-t(cbind(out$x,out$y))
}


#' @title Remove Alternate Values of a Matrix of Mid-infrared Spectra
#' @description Allow to remove alternate values of a matrix of Mid-Infrared Spectra
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param n numeric. Number of removed values between two variables. Defaults to 1.
#' @return A matrix spectra with removed values. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Removed spectra
#' rem <- sprem(spectra=spectra, n=1)
#' @export

sprem<-function(spectra, n=1){

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(n)) {
    stop('No n value provided')
  }

  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  k<-seq(from=1, to=dim(spectra)[2], by=n)
  output <-spectra[,k]
}


#' @title Normalizes the Absorbance Between 0 and 1 of a Matrix of Mid-infrared Spectra
#' @description Allows the normalization of the absorbance values between 0 and 1 of a matrix of mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @return A matrix spectra normalized between 0 and 1. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Normalized spectra between 0 and 1
#' norm <- spnorm01(spectra=spectra)
#' @export

spnorm01<-function(spectra){

  if (missing(spectra)) {
    stop('No spectral data provided')
  }

  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  s<-matrix(nrow=dim(spectra)[1]-1,ncol=dim(spectra)[2])
  for (i in 2:dim(spectra)[1]){
    s[i-1,]<-spectra[i,]-min(spectra[i,]) #offset
  }
  for (i in 2:dim(spectra)[1]){
    s[i-1,]<-spectra[i,]/max(spectra[i,])
  }
  output <- rbind(spectra[1, ],s)
}


#' @title Offset Correction of a Matrix of Mid-infrared Spectra
#' @description Allows the removal of the background of a matrix of mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @return A matrix spectra with with background values deleted. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Offset correction
#' offs <- spoffs(spectra=spectra)
#' @export

spoffs<-function(spectra){
  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  s<-matrix(nrow=dim(spectra)[1]-1,ncol=dim(spectra)[2])
  for (i in 2:dim(spectra)[1]){
    s[i-1,]<-spectra[i,]-min(spectra[i,])
  }
  output <- rbind(spectra[1, ],s)
}


#' @title Normalizes the Absorbance of a Matrix of Mid-infrared Spectra by a Specific Band
#' @description The absorbance values of the matrix of mid-infrared spectra is normalized by a specific band.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param wn numeric. Specific band (wavenumber) used to normalize the spectra.
#' @return A matrix spectra normalized by a specific band. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Normalized spectra by a specific band
#' normw <- spnormw(spectra=spectra, wn=1510)
#' @export

spnormw<-function(spectra, wn){
  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(wn)) {
    stop('No wl value provided')
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }

  index<-which(round(spectra[1, ])==wn, arr.ind=T)[1]
  s<-matrix(NA, nrow=nrow(spectra)-1,ncol=ncol(spectra))

  for(i in 1:nrow(s)) {
    s[i,]<-spectra[i+1,]/spectra[i+1,index]
  }

  output <- rbind(spectra[1, ],s)
}


#' @title Linear Baseline Correction of a Mid-infrared Spectrum
#' @description This function allows a linear correction of defects of the baseline of a mid-infrared spectrum.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectrum matrix. The matrix of FTIR spectrum. The first row corresponds to wavenumber; the second row corresponds to absorbance.
#' @param lbl vector. Vector of zero points of absorbance (two or more points).
#' @return A corrected spectrum matrix by means of a linear baseline. The first row corresponds to wavenumber; the second row corresponds to absorbance.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectrum <- as.matrix(t(spData[, c("Wavenumber","A")]))
#' # Linear baseline correction
#' mbl <- spmbl(spectrum=spectrum, lbl=c(1800, 1540, 840))
#' @export

spmbl<-function(spectrum,lbl){
  if (missing(spectrum)) {
    stop('No spectral data provided')
  }
  if (missing(lbl)) {
    stop('No vector of zero absorbance provided')
  }
  if (spectrum[1, 1] < spectrum[1, dim(spectrum)[2]]) {
    spectrum <- t(apply(spectrum, 1, rev))
  }
  if (lbl[1] < lbl[length(lbl)]) {
    lbl <- rev(lbl)
  }

  x<-spectrum[1,]
  y<-spectrum[2,]

  for (i in 1:length(lbl)){
    lbl[i]<-round(lbl[i], digits=0)
  }

  case <- 0
  if (which(round(spectrum[1,], digits=0)==lbl[1])[1]==1 & max(which(round(spectrum[1,], digits=0)==lbl[length(lbl)]))==length(x)){
    lbl <- lbl
    case <- 1
  }
  if (which(round(spectrum[1,], digits=0)==lbl[1])[1]==1 & max(which(round(spectrum[1,], digits=0)==lbl[length(lbl)]))!=length(x)){
    lbl <- c(lbl,round(x[length(x)], digits=0))
    case <- 2
  }
  if (which(round(spectrum[1,], digits=0)==lbl[1])[1]!=1 & max(which(round(spectrum[1,], digits=0)==lbl[length(lbl)]))==length(x)){
    lbl <- c(round(x[1], digits=0),lbl)
    case <- 3
  }
  if (which(round(spectrum[1,], digits=0)==lbl[1])[1]!=1 & max(which(round(spectrum[1,], digits=0)==lbl[length(lbl)]))!=length(x)){
    lbl <- c(round(x[1], digits=0),lbl,round(x[length(x)], digits=0))
    case <- 4
  }

  ap<- rep(0, length(lbl))
  for (i in 1:(length(lbl)-1)){
    if (i==1){
      ap[1] <- which(round(spectrum[1,], digits=0)==lbl[i], arr.ind=TRUE)[1]
      ap[length(lbl)] <-max(which(round(spectrum[1,], digits=0)==lbl[length(lbl)], arr.ind=TRUE))
    } else if((i>1)&(i<(length(lbl)))){
      ap[i] <-which(round(spectrum[1,], digits=0)==lbl[i], arr.ind=TRUE)[1]
    }
  }

  mp <- vector('list', length(ap)-1)
  if(length(ap)==2){
    ma<-matrix(NA, nrow=2, ncol=length(x))
    ma[1,]<-x[ap[1]:ap[2]]
    ma[2,]<-y[ap[1]:ap[2]]
    mp[[1]] <-ma
  }
  if(length(ap)>2){
    for (i in (1:(length(ap)-1))){
      if(i==1){
        ma<-matrix(NA, nrow=2, ncol=length(ap[1]:ap[2]))
        ma[1,]<-x[ap[1]:ap[2]]
        ma[2,]<-y[ap[1]:ap[2]]
        mp[[1]] <-ma
      }else {
        ma<-matrix(NA, nrow=2, ncol=length((ap[i]+1):ap[i+1]))
        ma[1,]<-x[(ap[i]+1):ap[i+1]]
        ma[2,]<-y[(ap[i]+1):ap[i+1]]
        mp[[i]]<-ma
      }
    }
  }

  if (case==1){
    m<- rep(NA, length(mp))
    for (i in 1:length(mp)){
      m[i]<-(mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
    }
    b<- rep(NA, length(mp))
    for (i in 1:length(mp)){
      b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
    }
  }
  if (case==2){
    m<- rep(NA, length(mp))
    for (i in 1:(length(mp))){
      if(i!=length(mp)){
        m[i]<-(mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
      }
      if(i==length(mp)){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          m[i]<- (mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
        }else{
          m[i]<-(mp[[i]][2,1]-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))])/(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))])
        }
      }
    }
    b<- rep(NA, length(mp))
    for (i in 1:length(mp)){
      if(i!=length(mp)){
        b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
      }
      if(i==length(mp)){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
        }else{
          b[i]<-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))]-m[i]*mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]
        }
      }
    }
  }

  if (case==3){
    m<- rep(NA, length(mp))
    for (i in 1:(length(mp))){
      if(i==1){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          m[i]<- (mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
        }else{
          m[i]<-(mp[[i]][2,1]-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))])/(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))])
        }
      }
      if(i!=1){
        m[i]<-(mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
      }
    }
    b<- rep(NA, length(mp))
    for (i in 1:length(mp)){
      if(i==1){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
        }else{
          b[i]<-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))]-m[i]*mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]
        }
      }
      if(i!=1){
        b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
      }
    }
  }


  if (case==4){
    m<- rep(NA, length(mp))
    for (i in 1:(length(mp))){
      if(i==1){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          m[i]<- (mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
        }else{
          m[i]<-(mp[[i]][2,1]-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))])/(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))])
        }
      }
      if(i>1&i<length(mp)){
        m[i]<-(mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
      }
      if(i==length(mp)){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          m[i]<-(mp[[i]][2,1]-mp[[i]][2,dim(mp[[i]])[2]])/(mp[[i]][1,1]-mp[[i]][1,dim(mp[[i]])[2]])
        }else{
          m[i]<-(mp[[i]][2,1]-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))])/(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))])
        }
      }
    }

    b<- rep(NA, length(mp))
    for (i in 1:length(mp)){
      if(i==1){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
        }else{
          b[i]<-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))]-m[i]*mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]
        }
      }
      if(i!=length(mp)){
        b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
      }
      if(i==length(mp)){
        if(mp[[i]][1,1]-mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]==0){
          b[i]<-mp[[i]][2,dim(mp[[i]])[2]]-m[i]*mp[[i]][1,dim(mp[[i]])[2]]
        }else{
          b[i]<-mp[[i]][2,which(mp[[i]][2,]==min(mp[[i]][2,]))]-m[i]*mp[[i]][1,which(mp[[i]][2,]==min(mp[[i]][2,]))]
        }
      }
    }
  }

  #Espectro corregido para cada intervalo
  lb <- vector('list', length(mp))
  for (i in 1:length(mp)){
    la<-matrix(NA, nrow=dim(mp[[i]])[1], ncol=dim(mp[[i]])[2])
    la[1,]<-mp[[i]][1,]
    for (j in 1:dim(mp[[i]])[2]){
      la[2,j]<-mp[[i]][2,j]-(mp[[i]][1,j]*m[i]+b[i])
    }
    lb[[i]]<-la
  }
  output <- as.matrix(data.frame(lb)) #joint list
}



#' @title Identification of Peaks of a Mid-infrared Spectra
#' @description This function allows to identify peaks of a mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param span numeric. Peak detection threshold.
#' @param tol numeric. Percentage of the maximum value of the spectrum (positive value).
#' @importFrom stats embed
#' @return An object of class sppeak, which is a list of matrices for each of the spectrum.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # List of peak detection
#' pks <- sppeak(spectra=spectra, span=3, tol=0.2)
#' # Peaks of the first spectrum
#' pks[[1]]
#' # Peaks of the second spectrum
#' pks[[2]]
#' @export

sppeak<-function(spectra, span=3, tol=0.2) {

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(span)) {
    stop('No span value provided')
  }
  if (missing(tol)) {
    stop('No tol value provided')
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }


  pklist<-list()
  for (i in 1:(nrow(spectra)-1)){
    span<-span
    tol<-tol
    dp<-as.numeric(spectra[i+1,])
    span.width <- span * 2 + 1
    loc.max <- span.width + 1 - apply(embed(dp, span.width), 1, which.max)
    loc.max[loc.max == 1 | loc.max == span.width] <- NA
    pk <- loc.max + 0:(length(loc.max) - 1)
    pk <- unique(pk[!is.na(pk)])

    nonda<-as.data.frame(spectra[1,])
    sp<-as.data.frame(spectra[i+1,])
    mPeak<-matrix(nrow = length(pk), ncol = 2)

    for (j in 1:length(pk)){
      mPeak[j,1]<-nonda[pk[j],1]
    }
    for (k in 1:length(pk)){
      mPeak[k,2]<-sp[pk[k],1]
    }
    if (tol>=0){
      mPeak<-mPeak[mPeak[,2]>max(mPeak[,2])*tol,]
    }else{
      mPeak<-mPeak[mPeak[,2]<min(mPeak[,2])*abs(tol),]
    }
    pklist[[i]] <- t(mPeak)
  }
  attr(pklist, "class") <- "peaks"
  pklist
}



#' @title Identification of Valleys of a Mid-infrared Spectra
#' @description This function allows to identify valleys of a mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param span numeric. Peak detection threshold.
#' @param tol numeric. Percentage of the maximum value of the spectrum (positive value).
#' @importFrom stats embed
#' @return An object of class spvalley, which is a list of matrices for each of the spectra.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # List of valley detection
#' vls <- sppeak(spectra=spectra, span=3, tol=0.2)
#' # Valleys of the first spectrum
#' vls[[1]]
#' # Valleys of the second spectrum
#' vls[[2]]
#' @export

spvalley<-function(spectra, span=3, tol=0.2){

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(span)) {
    stop('No span value provided')
  }
  if (missing(tol)) {
    stop('No tol value provided')
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }

  pklist<-list()
  tol=tol*(-1)
  for (i in 1:(nrow(spectra)-1))
  {

    dp<-as.numeric(spectra[i+1,])
    span.width <- span * 2 + 1
    loc.min <- span.width + 1 - apply(embed(dp, span.width), 1, which.min)
    loc.min[loc.min == 1 | loc.min == span.width] <- NA
    pk <- loc.min + 0:(length(loc.min) - 1)
    pk <- unique(pk[!is.na(pk)]) # PickPeaks

    nonda<-as.data.frame(spectra[1,])
    sp<-as.data.frame(spectra[i+1,])
    mPeak<-matrix(nrow = length(pk), ncol = 2) # definir matriz de peaks

    for (j in 1:length(pk)){
      mPeak[j,1]<-nonda[pk[j],1]
    }
    for (k in 1:length(pk)){
      mPeak[k,2]<-sp[pk[k],1]
    }
    if (tol>=0){
      mPeak<-mPeak[mPeak[,2]>max(mPeak[,2])*tol,]
    }else{
      mPeak<-mPeak[mPeak[,2]<min(mPeak[,2])*abs(tol),]
    }
    pklist[[i]] <- t(mPeak)
  }
  attr(pklist, "class") <- "peaks"
  pklist
}



#' @title Polynomial Baseline Correction of a Matrix of Mid-infrared Spectra
#' @description This function allows a polynomial correction of defects of the baseline of a mid-infrared spectrum.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param degree numeric. Degree of polynomial. Defaults to 2.
#' @param tol numeric. Tolerance of difference between iterations. Defaults to 0.001.
#' @param rep numeric. Maximum number of iterations. Defaults to 100.
#' @return An object of class spmblp, which is a list with the following components:
#'   \item{original}{Matrix of original mid-infrared spectra.}
#'   \item{baseline}{Matrix of polynomial baseline of mid-infrared spectra.}
#'   \item{corrected}{Matrix of polynomial baseline corrected mid-infrared spectra.}
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # List of polynomial baseline components
#' mblp <- spmblp(spectra, degree = 2, tol = 0.001, rep = 100)
#' # Original matrix
#' original <- mblp$original
#' # Baseline matrix
#' baseline <- mblp$baseline
#' # Corrected matrix
#' corrected <- mblp$corrected
#' @references
#' Lieber, C. A., and Mahadevan-Jansen, A. (2003). Automated method for subtraction of fluorescence from biological Raman spectra. Applied spectroscopy, 57(11), 1363-1367.
#' @export


spmblp <- function (spectra, degree = 2, tol = 0.001, rep = 100) {

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(degree)) {
    stop('No degree value provided')
  }
  if (missing(tol)) {
    stop('No tolerance value provided')
  }
  if (missing(rep)) {
    stop('No maximum number of iterations provided')
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }

  dimnames(spectra) <- NULL
  t<-spectra[1,]
  sp<-matrix(spectra[2:nrow(spectra),], nrow=nrow(spectra)-1, ncol=ncol(spectra))
  np <- dim(sp)

  baseline <- matrix(0, np[1], np[2])
  if (missing(t) || (t == FALSE))
    t <- 1:np[2]
  polx <- cbind(1/sqrt(np[2]), stats::poly(t, degree = degree))
  for (i in 1:np[1]) {
    ywork <- yold <- yorig <- sp[i, ]
    nrep <- 0
    repeat {
      nrep <- nrep + 1
      ypred <- polx %*% crossprod(polx, yold)
      ywork <- pmin(yorig, ypred)
      crit <- sum(abs((ywork - yold)/yold), na.rm = TRUE)
      if (crit < tol || nrep > rep)
        break
      yold <- ywork
    }
    baseline[i, ] <- ypred
  }
  value <-list(original = rbind(t, sp), baseline = rbind(t,baseline), corrected = rbind(t,sp - baseline))
  attr(value, "class") <- "spmblp"
  value
}



#' @title N-derived of a Mid-infrared Spectra
#' @description This function allows to determine the n-derivative of a mid-infrared spectra.
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param order numeric. Order of derivative. Defaults to 2.
#' @param p numeric. Polynomial order (p>order). Defaults to 3.
#' @param sw numeric. Filter length (must be odd). Defaults to 11.
#' @importFrom stats convolve
#' @importFrom pracma pinv
#' @return  A derivated spectra matrix. The first row corresponds to wavenumber; the second row corresponds to absorbance.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Derivative spectra
#' der <- spder(spectra=spectra, order=2, p=3, sw= 11)
#' @export


spder<-function (spectra, order=2, p=3, sw= 11) {

  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(order)) {
    stop('No order of derivate provided')
  }
  if (missing(p)) {
    stop('No polynomial order provided')
  }
  if (missing(sw)) {
    stop('No filter length provided')
  }
  stopifnot(is.numeric(spectra), is.numeric(sw))
  if (sw <= 1 || sw%%2 == 0) {
    stop("Argument sw must be an odd integer greater than 1.")
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }
  st<-matrix(NA, nrow=nrow(spectra)-1, ncol=ncol(spectra))
  for (i in 1:(nrow(spectra)-1)){
    # Inserting on both sides N points with intensity equal to the first and last value
    s <- spectra[(i+1),] #Data
    sa <- rep(NA, length(s)+2*sw)
    s <- c(rep(s[1],sw),s,rep(s[length(s)],sw))

    fc <- (sw - 1)/2
    X <- outer(-fc:fc, 0:p, FUN = "^")
    Y <- pinv(X)
    t <- convolve(s, rev(Y[(order + 1), ]), type = "o")
    Ts <- t[(fc + 1):(length(t) - fc)]
    st[i,]<-Ts[(sw+1):(length(Ts)-sw)]
  }
  output <- rbind(spectra[1,], st)
}


#' @title Savitzky-Golay Smoothing Filter of a Mid-infrared Spectrum
#' @description This function allows applying a Savitzky-Golay smoothing filter to the mid-infrared spectrum (N spectra= 1).
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectrum matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the second row corresponds to absorbances.
#' @param p numeric. Filter order. Defaults to 2.
#' @param sw numeric. Filter length (must be odd). Defaults to 21.
#' @importFrom stats convolve
#' @importFrom pracma pinv
#' @return A smoothed spectrum matrix by means of a Savitzky-Golay smoothing filter. The first row corresponds to wavenumber; the second row corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectrum <- as.matrix(t(spData[, c("Wavenumber","A")]))
#' # Smoothed spectrum
#' sga <- spsga(spectrum=spectrum, p=2, sw= 21)
#' @export


spsga<-function (spectrum, p=2, sw= 21) {

  if (missing(spectrum)) {
    stop('No spectruml data provided')
  }
  if (missing(sw)) {
    stop('No filter length provided')
  }
  if (missing(p)) {
    stop('No filter order provided')
  }
  stopifnot(is.numeric(spectrum), is.numeric(sw))
  if (sw <= 1 || sw%%2 == 0) {
    stop("Argument 'sw' must be an odd integer greater than 1.")
  }
  if (spectrum[1, 1] < spectrum[1, dim(spectrum)[2]]) {
    spectrum <- t(apply(spectrum, 1, rev))
  }

  # Inserting on both sides N points
  s <- spectrum[2,] #Data
  sa <- rep(NA, length(s)+2*sw)
  s <- c(rep(s[1],sw),s,rep(s[length(s)],sw))

  fc <- (sw - 1)/2
  X  <- outer(-fc:fc, 0:p, FUN = "^")
  Y  <- pinv(X)
  t  <- convolve(s, rev(Y[1,]), type = "o")
  Ts <- t[(fc + 1):(length(t) - fc)]

  output <- rbind(spectrum[1,], Ts[(sw+1):(length(Ts)-sw)])
}


#' @title Savitzky-Golay Smoothing Filter of a Mid-infrared Spectra
#' @description This function allows applying a Savitzky-Golay smoothing filter to the mid-infrared spectra (N spectra > 1).
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectra matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @param p numeric. Filter order. Defaults to 2.
#' @param sw numeric. Filter length (must be odd). Defaults to 21.
#' @importFrom stats convolve
#' @importFrom pracma pinv
#' @return A smoothed spectra matrix by means of a Savitzky-Golay smoothing filter. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectra <- as.matrix(t(spData))
#' # Smoothed spectra
#' sgb <- spsgb(spectra=spectra, p=2, sw= 21)
#' @export


spsgb<-function (spectra, p=2, sw= 21) {


  if (missing(spectra)) {
    stop('No spectral data provided')
  }
  if (missing(p)) {
    stop('No polynomial order provided')
  }
  if (missing(sw)) {
    stop('No filter length provided')
  }
  stopifnot(is.numeric(spectra), is.numeric(sw))
  if (sw <= 1 || sw%%2 == 0) {
    stop("Argument 'sw' must be an odd integer greater than 1.")
  }
  if (spectra[1, 1] < spectra[1, dim(spectra)[2]]) {
    spectra <- t(apply(spectra, 1, rev))
  }

  st<-matrix(NA, nrow=nrow(spectra)-1, ncol=ncol(spectra))
  for (i in 1:(nrow(spectra)-1)){
    # Inserting on both sides N points with intensity equal to the first and last value
    s <- spectra[(i+1),] #Data
    sa <- rep(NA, length(s)+2*sw)
    s <- c(rep(s[1],sw),s,rep(s[length(s)],sw))

    fc <- (sw - 1)/2
    X  <- outer(-fc:fc, 0:p, FUN = "^")
    Y  <- pinv(X)
    t  <- convolve(s, rev(Y[1,]), type = "o")
    Ts <- t[(fc + 1):(length(t) - fc)]
    st[i,]<-Ts[(sw+1):(length(Ts)-sw)]
  }
  output <- rbind(spectra[1,], st)
}



#' @title  Moving-average Smoothing Filter of a Mid-infrared Spectrum
#' @description This function allows applying a Moving-average smoothing filter to the mid-infrared spectrum (N spectra = 1).
#' @author Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
#' @param spectrum matrix. The matrix of FTIR spectrum. The first row corresponds to wavenumber; the second row corresponds to absorbances.
#' @param N numeric. Length of the smoothing window. Defaults to 21.
#' @importFrom stats convolve
#' @importFrom pracma pinv
#' @return A smoothed spectrum matrix by means of a Moving-average smoothing filter. The first row corresponds to wavenumber; the second row corresponds to absorbances.
#' @examples
#' data(spData)
#' # Convert data frame to matrix
#' spectrum <- as.matrix(t(spData[, c("Wavenumber","A")]))
#' # Smoothed spectrum
#' mws <- spmws(spectrum = spectrum, N = 21)
#' @export


spmws<-function (spectrum, N=21) {

  if (missing(spectrum)) {
    stop('No spectral data provided')
  }
  if (missing(N)) {
    stop('No length of the smoothing window provided')
  }
  stopifnot(is.numeric(spectrum), is.numeric(N))
  if (spectrum[1, 1] < spectrum[1, dim(spectrum)[2]]) {
    spectrum <- t(apply(spectrum, 1, rev))
  }

  # Inserting on both sides N points
  s <- spectrum[2,] #Data
  sa <- rep(NA, length(s)+2*N)
  sa <- c(rep(s[1],N),s,rep(s[length(s)],N))

  # Replacement of the kth point by the average of its 2N next neighbours (N to the left, N to the right)
  for(k in N:(length(sa)-N)){
    sa[k]<-(1/(2*N))*(sum(sa[(k-N):(k-1)])+sum(sa[(k+1):(k+N)]))
  }

  output <- as.matrix(rbind(spectrum[1,], sa[(N+1):(length(sa)-N)]))
}


