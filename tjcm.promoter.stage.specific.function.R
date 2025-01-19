REP = function (D, k) 
{
  r <- nrow(D)
  c <- ncol(D)
  DD <- NULL
  for (i in 1:c) {
    DDi <- rep(D[, i], k)
    DD <- cbind(DD, DDi)
  }
  colnames(DD) <- colnames(D)
  as.data.frame(DD)
}
Formula0 = function (names) 
{
  formula <- "y~"
  if (length(names) == 1) {
    formula = paste(formula, names[1], "+ transcript")
  }
  else if (length(names) > 1) {
    for (i in 1:(length(names))) {
      formula <- paste(formula, names[i], "+")
    }
    formula <- paste(formula, "transcript")
  }
  formula <- as.formula(formula)
  formula
}

Formula1 = function (names) 
{
  formula <- "y~"
  if (length(names) == 1) {
    formula = paste(formula, names[1], "* transcript")
  }
  else if (length(names) > 1) {
    formula <- paste(formula, "(")
    for (i in 1:(length(names) - 1)) {
      formula <- paste(formula, names[i], "+")
    }
    formula <- paste(formula, names[length(names)])
    formula <- paste(formula, ") * transcript")
  }
  formula <- as.formula(formula)
  formula
}

negative.binomial = function (theta = stop("'theta' must be specified"), link = "log") 
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt")) {
    stats <- make.link(linktemp)
  } else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    } else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"", 
                         linktemp))
  }
  .Theta <- theta
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", theta, envir = env)
  variance <- function(mu) mu + mu^2/.Theta
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) + 
      lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - 
      lgamma(.Theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rnegbin(nsim * length(ftd), ftd, .Theta)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("Negative Binomial(", format(round(theta, 
                                                      4)), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun), 
            class = "family")
}

p.vector = function (data, design, Q = 0.05, MT.adjust = "BH", min.obs = 6,
                     counts = FALSE, family = NULL, theta = 10, epsilon = 1e-05,
                     item = "gene")
{
  if (is.data.frame(design) || is.matrix(design)) {
    dis <- design
    groups.vector = NULL
    edesign = NULL
  } else if (is.list(design)) {
    dis <- as.data.frame(design$dis)
    groups.vector <- design$groups.vector
    edesign <- design$edesign
  }
  if (is.null(family)) {
    if (!counts) {
      family = gaussian()
    }
    if (counts) {
      family = negative.binomial(theta)
    }
  }
  dat <- as.matrix(data)
  dat <- dat[, as.character(rownames(dis))]
  G <- nrow(dat)
  count.na <- function(x) (length(x) - length(x[is.na(x)]))
  dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
  sumatot <- apply(dat, 1, sum)
  counts0 <- which(sumatot == 0)
  if (length(counts0) > 0) {dat <- dat[-counts0, ]}
  g <- dim(dat)[1]
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  p.vector <- vector(mode = "numeric", length = g)
  for (i in 1:g) {
    y <- as.numeric(dat[i, ])
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) print(paste(c("fitting ", item, i, "out of", g), collapse = " "))
    tryCatch(
      {
        model.glm <- glm(y ~ ., data = dis, family = family, epsilon = epsilon) ##calculate coefficient for time1,2,3, fit a linear model
        if (model.glm$null.deviance == 0) {
          p.vector[i] = 1
        } else {
          model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon) ##fit the model to the mean value
          if (family$family == "gaussian") {
            test <- anova(model.glm.0, model.glm, test = "F")
            if (is.na(test[6][2, 1])) {
              p.vector[i] = 1
            } else {
              p.vector[i] = test[6][2, 1]
            }
          } else {
            test <- anova(model.glm.0, model.glm, test = "Chisq")
            if (is.na(test[5][2, 1])) {
              p.vector[i] = 1
            } else {
              p.vector[i] = test[5][2, 1]
            }
          }
        }
      },error=function(cond){p.vector[i] = 1}
    )

  }
  p.adjusted <- p.adjust(p.vector, method = MT.adjust, n = length(p.vector))
  genes.selected <- rownames(dat)[which(p.adjusted <= Q)]
  FDR <- sort(p.vector)[length(genes.selected)]
  SELEC <- as.matrix(as.data.frame(dat)[genes.selected, ])
  if (nrow(SELEC) == 0) print("no significant genes")
  p.vector <- as.matrix(p.vector)
  rownames(p.vector) <- rownames(dat)
  colnames(p.vector) <- c("p.value")
  output <- list(SELEC, p.vector, p.adjusted, G, g, FDR, nrow(SELEC),
                 dis, dat, min.obs, Q, groups.vector, edesign, family)
  names(output) <- c("SELEC", "p.vector", "p.adjusted", "G",
                     "g", "FDR", "i", "dis", "dat", "min.obs", "Q", "groups.vector",
                     "edesign", "family")
  output
}

T.fit = function (data, design = data$dis, step.method = "backward", 
                  min.obs = data$min.obs, alfa = data$Q, nvar.correction = FALSE, 
                  family = gaussian(), epsilon = 1e-05, item = "gene") 
{
  if (is.list(data)) {
    dat <- as.matrix(data$SELEC)
    dat <- rbind(c(rep(1, ncol(dat))), dat)
    groups.vector <- data$groups.vector
    groups.vector <- c(groups.vector[nchar(groups.vector) == 
                                       min(nchar(groups.vector))][1], groups.vector)
    edesign <- data$edesign
    G <- data$g
    family <- data$family
  } else {
    G <- nrow(data)
    data <- rbind(c(rep(1, ncol(data))), data)
    dat <- as.matrix(data)
    count.na <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[apply(dat, 1, count.na) >= min.obs, ]
    groups.vector = NULL
    edesign = NULL
  }
  dis <- as.data.frame(design)
  dat <- dat[, as.character(rownames(dis))]
  g <- (dim(dat)[1] - 1)
  n <- dim(dat)[2]
  p <- dim(dis)[2]
  vars.in <- colnames(dis)
  sol <- coefficients <- group.coeffs <- t.score <- sig.profiles <- NULL
  influ.info <- matrix(NA, nrow = nrow(dis), ncol = 1)
  rownames(influ.info) <- rownames(dis)
  if (nvar.correction) alfa <- alfa/ncol(dis)
  for (i in 2:(g + 1)) {
    y <- as.numeric(dat[i, ])
    name <- rownames(dat)[i]
    reg <- try(stepback(y = y, d = dis, alfa = alfa, family = family, epsilon = epsilon))
    if (class(reg) == "try-error") {
      reg <- try(stepfor(y = y, d = dis, alfa = alfa, family = family, epsilon = epsilon))
    }
    if (class(reg) == "try-error") {
      reg <- try(two.ways.stepback(y = y, d = dis, alfa = alfa,family = family, epsilon = epsilon))
    }
    if (class(reg) == "try-error") {
      reg <- try(two.ways.stepfor(y = y, d = dis, alfa = alfa,family = family, epsilon = epsilon))
    }
    if (class(reg) == "try-error") {
      reg <- try(stepback(y = y, d = dis, alfa = alfa, family = gaussian(), epsilon = epsilon))
    }
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) print(paste(c("fitting ", item, i, "out of", g),collapse = " "))
    lmf <- try(glm(y ~ ., data = as.data.frame(dis), family = family, epsilon = epsilon))
    if (class(lmf) == "try-error") {lmf <- glm(y ~ ., data = as.data.frame(dis), family = gaussian(), epsilon = epsilon)}
    result <- summary(lmf)
    novar <- vars.in[!is.element(vars.in, names(result$coefficients[, 4]))]
    influ <- influence.measures(reg)$is.inf
    influ <- influ[, c(ncol(influ) - 3, ncol(influ) - 1)]
    influ1 <- which(apply(influ, 1, all))
    if (length(influ1) != 0) {
      paste.names <- function(a) {
        paste(names(a)[a], collapse = "/")
      }
      match <- match(rownames(dis), rownames(influ))
      influ <- as.data.frame(apply(influ, 1, paste.names))
      influ.info <- cbind(influ.info, influ[match, ])
      colnames(influ.info)[ncol(influ.info)] <- name
    }
    result <- summary(reg)
    if ((!(result$aic == -Inf) & !is.na(result$aic) & family$family == 
         "gaussian") | family$family != "gaussian") {
      k <- i
      model.glm.0 <- glm(y ~ 1, family = family, epsilon = epsilon)
      if (family$family == "gaussian") {
        test <- anova(model.glm.0, reg, test = "F")
        p.value = test[6][2, 1]
      } else {
        test <- anova(model.glm.0, reg, test = "Chisq")
        p.value = test[5][2, 1]
      }
      bondad <- (reg$null.deviance - reg$deviance)/reg$null.deviance
      if (bondad < 0) {
        bondad = 0
      }
      beta.coeff <- result$coefficients[, 1]
      beta.p.valor <- result$coefficients[, 4]
      coeff <- rep(0, (length(vars.in) + 1))
      if (length(novar) != 0) {
        for (m in 1:length(novar)) {
          coeff[position(dis, novar[m]) + 1] <- NA
        }
      }
      p.valor <- t <- as.numeric(rep(NA, (length(vars.in) + 1)))
      if (result$coefficients[, 4][rownames(result$coefficients) == "(Intercept)"] < alfa) {
        coeff[1] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                               "(Intercept)"]
        p.valor[1] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                 "(Intercept)"]
        t[1] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                           "(Intercept)"]
      }
      for (j in 2:length(coeff)) {
        if (is.element(vars.in[j - 1], rownames(result$coefficients))) {
          coeff[j] <- result$coefficients[, 1][rownames(result$coefficients) == 
                                                 vars.in[j - 1]]
          p.valor[j] <- result$coefficients[, 4][rownames(result$coefficients) == 
                                                   vars.in[j - 1]]
          t[j] <- result$coefficients[, 3][rownames(result$coefficients) == 
                                             vars.in[j - 1]]
        }
      }
      if (!all(is.na(p.valor))) {
        sol <- rbind(sol, as.numeric(c(p.value, bondad, 
                                       p.valor)))
        coefficients <- rbind(coefficients, coeff)
        t.score <- rbind(t.score, t)
        sig.profiles <- rbind(sig.profiles, y)
        h <- nrow(sol)
        rownames(sol)[h] <- name
        rownames(coefficients)[h] <- name
        rownames(t.score)[h] <- name
        rownames(sig.profiles)[h] <- name
      }
    }
  }
  message("step wise regression finished")
  if (!is.null(sol)) {
    sol <- as.data.frame(sol)
    coefficients <- as.data.frame(coefficients)
    coeffic <- coefficients
    t.score <- as.data.frame(t.score)
    sig.profiles <- as.data.frame(sig.profiles)
    colnames(sol) <- c("p-value", "R-squared", "p.valor_beta0", 
                       paste("p.valor_", vars.in, sep = ""))
    colnames(coefficients) <- c("beta0", paste("beta", vars.in, 
                                               sep = ""))
    colnames(t.score) <- c("t.score_beta0", paste("t.score_", 
                                                  vars.in, sep = ""))
    colnames(sig.profiles) <- colnames(dat)
    if (!is.null(groups.vector) & !is.null(edesign)) {
      groups <- colnames(edesign)[3:ncol(edesign)]
      degree <- (length(groups.vector)/length(groups)) - 1
      for (w in 1:nrow(coefficients)) {
        A <- NULL
        col.names <- NULL
        for (l in 1:length(groups)) {
          B <- reg.coeffs(coefficients = coefficients[w, ], groups.vector = groups.vector, group = groups[l])
          cols <- paste(rep(groups[l], each = length(B)), paste("beta", c(0:(length(B) - 1)), sep = ""), sep = "_")
          A <- c(A, B)
          col.names <- c(col.names, cols)
        }
        group.coeffs <- (rbind(group.coeffs, A))
      }
      colnames(group.coeffs) <- col.names
      rownames(group.coeffs) <- rownames(coefficients)
    }
  }
  if (ncol(influ.info) > 2) {
    print(paste("Influence:", ncol(influ.info) - 1, "genes with influential data at slot influ.info. Model validation for these genes is recommended"))
  }
  influ.info <- influ.info[, -1]
  output <- list(sol, sig.profiles, coefficients, as.data.frame(group.coeffs), 
                 t.score, vars.in, G, g, dat, dis, step.method, groups.vector, 
                 edesign, influ.info)
  names(output) <- c("sol", "sig.profiles", "coefficients", 
                     "group.coeffs", "t.score", "variables", "G", "g", "dat", 
                     "dis", "step.method", "groups.vector", "edesign", "influ.info")
  output
}


IsoModel = function(data, gen, design = NULL, Q = 0.05, Qfit = 0.05, min.obs = 6, minorFoldfilter = NULL, 
                    counts = FALSE, family = NULL, theta = 10, epsilon = 1e-05, step.method = "forward") 
{
  Genes <- unique(gen)
  g <- length(Genes)
  if (is.null(family)) {
    if (!counts) {
      family = gaussian()
    }
    if (counts) {
      family = negative.binomial(theta)
    }
  }
  print(paste(nrow(data), "transcripts"))
  print(paste(length(unique(gen)), "genes"))
  if (!is.null(minorFoldfilter)) {
    print("Removing low expressed minor isoforms")
    moreOne <- names(which(table(gen) > 1))
    iso.sel <- NULL
    gene.sel <- NULL
    for (i in moreOne) {
      which(gen == i)
      gene.data <- data[which(gen == i), ]
      isoSUM <- apply(gene.data, 1, sum)
      major <- names(which(isoSUM == max(isoSUM)))[1]
      minors <- names(which(isoSUM != max(isoSUM)))
      div <- as.numeric(matrix(rep(gene.data[major, ],length(minors)), ncol = ncol(data), length(minors), 
                               byrow = T))/as.matrix(gene.data[minors, ])
      is <- names(which(apply(div, 1, min, na.rm = T) < minorFoldfilter))
      iso.sel <- c(iso.sel, is)
      gene.sel <- c(gene.sel, rep(i, length(is)))
    }
    data <- data[iso.sel, ]
    gen <- gene.sel
    print(dim(data))
    print(length(gen))
    print("Done")
    print(paste(nrow(data), "remaining transcripts"))
    print(paste(length(unique(gen)), "remaining genes"))
  }
  NT <- tapply(rownames(data), gen, length)
  Genes1 <- names(which(NT != 1))
  data1 <- data[gen %in% Genes1, ]
  gen1 <- gen[gen %in% Genes1]
  Genes1 <- unique(gen1)
  g <- length(Genes1)
  print(paste("There are",g, "genes have >=2 transcripts"))
  dis <- as.data.frame(design$dis)
  mycolnames <- colnames(dis)
  pval <- NULL
  for (i in c(1:g)) {
    div <- c(1:round(g/100)) * 100
    if (is.element(i, div)) print(paste(c("fitting gene", i, "out of", g), collapse = " "))
    zz <- data1[gen1 == Genes1[i], ]
    nt <- nrow(zz)
    dis.gen <- REP(dis, nt)
    y <- c(t(as.matrix(zz)))
    transcript <- factor(rep(c(1:nt), each = ncol(zz)))
    ydis <- cbind(y, dis.gen, transcript)
    pvali = tryCatch(
      {
        model0 <- glm(Formula0(mycolnames), data = ydis, family = family, epsilon = epsilon)
        model1 <- glm(Formula1(mycolnames), data = ydis, family = family, epsilon = epsilon)
        if (family$family == "gaussian") {
           pvali = anova(model0, model1, test = "F")[2, 6]
        } else {
          pvali = anova(model0, model1, test = "Chisq")[2,5]
        }
      },
      error=function(cond){return(1)}
    )
    names(pvali) <- Genes1[i]
    pval <- c(pval, pvali)
  }
  num.genes <- sum(p.adjust(pval) <= Q, na.rm = TRUE) ##only keep genes that are significantly expressed
  selected.genes <- names(sort(p.adjust(pval))[1:num.genes])
  data2 <- data[gen %in% selected.genes, ]
  gen2 <- gen[gen %in% selected.genes]
  ## calculate p-value to Finding significant isoforms
  pvector2 <- p.vector(data=data2, design=design, counts = counts, item = "isoform", Q = Qfit, theta=theta)
  ## calculate R correlation to filter significant genes, and if there are several groups, here is to find differences between experimental groups
  message("starting Tfit")
  Tfit2 <-T.fit(data = pvector2, step.method = step.method, item = "isoform") #, family = pvector2$family
  ISO.SOL <- list(data, gen, design, selected.genes, pvector2, Tfit2)
  names(ISO.SOL) <- c("data", "gen", "design", "DSG", "pvector", "Tfit")
  ISO.SOL
}

getDS = function (Model, vars = "all", rsq = 0.4) 
{
  data <- Model$data
  gen <- Model$gen
  selected.genes <- Model$DSG
  Tfit2 <- Model$Tfit
  data2 <- data[gen %in% selected.genes, ]
  gen2 <- gen[gen %in% selected.genes]
  get2 <- get.siggenes(tstep = Tfit2, vars = vars, rsq = rsq)
  sig.iso2 <- get2$summary
  if(!is.null(dim(sig.iso2)) & vars=="groups") {
    sig.iso2 = as.character(get2$summary[,1])
  } else if (!is.null(dim(sig.iso2)) & vars=="each") {
    sig.iso2 = as.character(get2$summary[,-1])}
  gen.sig.iso2 <- as.character(gen2[rownames(data2) %in% sig.iso2])
  NumIso.by.gene <- tapply(sig.iso2, gen.sig.iso2, length)
  DSG_distributed_by_number_of_DETs <- NumIso.by.gene
  T.iso2 <- table(DSG_distributed_by_number_of_DETs, useNA = "ifany")
  print(paste(length(selected.genes), " DSG selected"))
  print(paste(length(gen.sig.iso2), " DETs selected"))
  print("Statistics of each gene's DET number")
  print(T.iso2)
  List0 <- setdiff(selected.genes, gen.sig.iso2)
  out <- list(Model, get2, selected.genes, sig.iso2, List0,NumIso.by.gene)
  names(out) <- c("Model", "get2", "DSG", "DET", "List0", "NumIso.by.gene")
  out
}

seeDS = function (get, rsq = 0.4, cluster.all = TRUE, plot.mDSG = FALSE, 
                  k = 6, cluster.method = "hclust", k.mclust = FALSE) 
{
  Model <- get$Model
  get2 <- get$get2
  data <- Model$data
  gen <- Model$gen
  design <- Model$design
  sig.iso2 <- get2$summary
  gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
  NT2 <- get$NumIso.by.gene
  if (cluster.all) {
    step3 <- p.vector(data, design, family = Model$pvector2$family)
    Tfit3 <- T.fit(step3)
    get3 <- get.siggenes(Tfit3, vars = "all", rsq = rsq)
    sig.iso3 <- get3$summary
    H <- see.genes(data = get3$sig.genes, item = "Isoforms", cluster.method = cluster.method, k.mclust=k.mclust, k = k, newX11 = FALSE)
    cut <- H$cut[sig.iso2]
  }
  if (!cluster.all) {
    H <- see.genes(data = get2$sig.genes, item = "Isoforms", cluster.method = cluster.method, k.mclust=k.mclust, k = k, newX11 = FALSE)
    cut <- H$cut
  }
  if (plot.mDSG) {
    data.clust <- get2$sig.genes$sig.profiles
    genes.1 <- names(NT2[NT2 == 1])
    data.clust1 <- data.clust[gen.sig.iso2 %in% genes.1,]
    H1 <- see.genes(data = data.clust1, edesign = design$edesign, 
                    cluster.method = cluster.method, k.mclust = k.mclust, 
                    k = k, item = "Isoforms", newX11 = FALSE)
  }
  out <- list(Model, get2, NT2, cut, gen.sig.iso2)
  names(out) <- c("Model", "get2", "NumIso.by.gene", "cut", 
                  "names.genes")
  out
}

tableDS = function (seeDS)
{
  Model <- seeDS$Model
  get2 <- seeDS$get2
  cut <- seeDS$cut
  data <- Model$data
  gen <- Model$gen
  design <- Model$design
  sig.iso2 <- get2$summary
  gen.sig.iso2 <- seeDS$names.genes #as.character(gen[rownames(data) %in% sig.iso2])
  NT2 <- seeDS$NumIso.by.gene
  data.clust <- get2$sig.genes$sig.profiles
  genes.2 <- names(NT2[NT2 > 1])
  data.clust <- data.clust[gen.sig.iso2 %in% genes.2, ]
  gen.sig.iso22 <- gen.sig.iso2[gen.sig.iso2 %in% genes.2]
  cut <- cut[sig.iso2[gen.sig.iso2 %in% genes.2]]
  gen.sig.iso2 <- gen.sig.iso22
  unic <- unique(gen.sig.iso2)
  G2 <- length(unic)
  Mayor = NULL
  for (j in 1:G2) {
    zz <- data.clust[gen.sig.iso2 == unic[j], ]
    Mayor.j <- MayorIso(zz)
    names(Mayor.j) <- rownames(zz)
    Mayor <- c(Mayor, Mayor.j)
  }
  Mayor = Mayor[names(cut)]
  cuts <- NULL
  cuts.name = NULL
  for (i in 1:length(unic)) {
    cutMi <- cut[gen.sig.iso2 == unic[i] & Mayor == 1]
    cutmi <- sort(unique(cut[gen.sig.iso2 == unic[i] & Mayor == 0]))
    cutMi.name = names(cutMi)
    cutmi.name = names(cut[gen.sig.iso2 == unic[i] & Mayor == 0])
    if (length(cutmi.name) > 1) {
      cutmi <- paste(cutmi, collapse = "_")
      cutmi.name = paste(cutmi.name, collapse = "_")
    }
    cuti <- c(cutMi, cutmi)
    cuts <- rbind(cuts, cuti)
    rownames(cuts)[i] <- unic[i]
    
    cuti.name <- c(cutMi.name, cutmi.name)
    cuts.name <- rbind(cuts.name, cuti.name)
    rownames(cuts.name)[i] <- unic[i]
    
  }
  cuts <- as.data.frame(cuts)
  colnames(cuts) <- c("Cluster.Mayor", "Cluster.minor")
  cuts.name <- as.data.frame(cuts.name)
  colnames(cuts.name) = c("Major", "Minor")
  IsoTable <- table(cuts)
  out <- list(IsoTable, cuts, cuts.name)
  names(out) <- c("IsoTable", "IsoClusters", "IsoMajorMinor")
  out
}

MayorIso = function (zz) 
{
  if (is.null(nrow(zz))) { sol = 1 }
  else { M <- apply(zz, 1, sum)
  sol = as.numeric(M == max(M)) }
  sol
}

PodiumChange = function (get, only.sig.iso = FALSE, comparison = c("any", "groups","specific"), group.name = "Ctr", time.points = 0) 
{
  Model <- get$Model
  get2 <- get$get2
  data <- Model$data
  gen <- Model$gen
  edesign <- Model$design$edesign
  repvect = edesign[, 2]
  if (only.sig.iso) {
    sig.iso2 <- get2$summary
    gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
    NT2 <- get$NumIso.by.gene
    data.clust <- as.matrix(get2$sig.genes$sig.profiles)
    genes.2 <- names(NT2[NT2 > 1])
    data.clust <- data.clust[gen.sig.iso2 %in% genes.2, ]
    gen.sig.iso22 <- gen.sig.iso2[gen.sig.iso2 %in% genes.2]
    gen.sig.iso2 <- gen.sig.iso22
  } else {
    sig.iso2 <- get2$summary
    gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
    data.clust <- as.matrix(data[gen %in% unique(gen.sig.iso2), ]) #Model$DSG
    sig.iso2 <- rownames(data.clust)
    gen.sig.iso2 <- as.character(gen[rownames(data) %in% sig.iso2])
    NT2 <- tapply(sig.iso2, gen.sig.iso2, length)
  }
  time.M <- tapply(edesign[, 1], repvect, mean)
  groups.M <- apply(edesign[, 3:ncol(edesign),drop=FALSE], 2, function(x) { tapply(x, repvect, mean)})
  unic <- unique(gen.sig.iso2)
  Mayor = NULL
  LIST = NULL
  for (i in 1:length(unic)) {
    zz <- data.clust[gen.sig.iso2 == unic[i], ]
    M <- MayorIso(zz)
    if(sum(M)>1){M[M==1][-1] = 0} ##if two isoform is exactly the same expression, use the first one as the major isoform
    zzM <- zz[M == 1, ]
    MzzM <- tapply(zzM, repvect, mean)
    zzm <- zz[M != max(M), ]
    if (is.null(nrow(zzm))) ni = 1 else ni = nrow(zzm)
    if (ni == 1) {Mzzm = tapply(zzm, repvect, mean)
    } else {Mzzm <- t(apply(zzm, 1, function(x) { tapply(x, repvect, mean) }))}
    if (ni == 1) dif = MzzM - Mzzm else dif <- t(apply(Mzzm, 1, function(x) { MzzM - x }))
    if (comparison == "any") {
      if (any(dif < 0)) LIST <- c(LIST, unic[i])
    } else if (comparison == "specific") {
      col <- groups.M[, colnames(groups.M) == group.name]
      if (ni == 1) {change <- all(dif[col == 1 & time.M == time.points] < 0)
      } else change <- apply(dif[, col == 1 & time.M == time.points], 1, function(x) { all(x < 0) })
      if (any(change))  LIST <- c(LIST, unic[i])
    } else if (comparison == "group") {
      mayors = NULL
      for (k in 3:ncol(edesign)) {
        mayors = cbind(mayors, MayorIso(zz[, edesign[, k] == 1]))
      }
      if (all(mayors - mayors[, 1] != 0)) 
        LIST <- c(LIST, unic[i])
    }
  }
  gen.L <- gen.sig.iso2[gen.sig.iso2 %in% LIST]
  data.L <- data.clust[gen.sig.iso2 %in% LIST, ]
  output <- list(LIST, data.L, gen.L, edesign)
  names(output) <- c("L", "data.L", "gen.L", "edesign")
  output
}

PlotGroups = function (data, edesign = NULL, time = edesign[, 1], groups = edesign[, c(3:ncol(edesign))], 
                       repvect = edesign[, 2], show.fit = FALSE, x.labels = NULL,
                       dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, 
                       nvar.correction = FALSE, summary.mode = "median", show.lines = TRUE, 
                       groups.vector = NULL, xlab = "Time", ylab = "Expression value", 
                       cex.xaxis = 1, ylim = NULL, main = NULL, cexlab = 0.8, legend = TRUE, 
                       sub = NULL, item = NULL) 
{
  time.change = 1:length(unique(time))
  names(time.change) = unique(time)
  time = as.numeric(time.change[as.character(time)])
  if (!is.vector(data)) {
    if (summary.mode == "representative") {
      distances <- apply(as.matrix(dist(data, diag = TRUE, upper = TRUE)), 1, sum)
      representative <- names(distances)[distances == min(distances)]
      yy <- as.numeric(data[rownames(data) == representative, ])
      sub <- paste("Representative:", representative)
    } else if (summary.mode == "median") {
      yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Median profile of", nrow(data), item, sep = " ")
      }
    } else if (summary.mode == "mean") {
      yy <- apply(as.matrix(data), 2, mean, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Mean profile of", nrow(data), item, sep = " ")
      }
    } else stop("not valid summary.mode")
    if (dim(data)[1] == 1) {
      sub <- rownames(data)
    }
  } else if (length(data) != 0) {
    yy <- as.numeric(data)
    sub <- rownames(data)
  } else stop("empty data")
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  reps <- i.rank(repvect)
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y)) {
    y[k] <- mean(yy[reps == k], na.rm = TRUE)
    x[k] <- mean(time[reps == k])
    for (j in 1:ncol) {
      g[k, j] <- mean(groups[reps == k, j])
    }
  }
  if (is.null(ylim)) ylim = c(min(as.numeric(yy), na.rm = TRUE), max(as.numeric(yy), na.rm = TRUE))
  abcissa <- x
  xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE) * 1.3)
  color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
  color2 <- groups
  for (j in 1:ncol) {
    color2[, j] <- color2[, j] * j
  }
  color2 <- as.vector(apply(color2, 1, sum) + 1)
  plot(x = time, y = yy, pch = 21, xlab = xlab, ylab = ylab, 
       xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim, 
       cex = cexlab, col = color2)
  if (is.null(x.labels)){
    axis(1, at = unique(abcissa), labels = x.labels, cex.axis = cex.xaxis)
  } else {
    axis(1, at = unique(abcissa), labels = unique(abcissa), cex.axis = cex.xaxis)
  }
  if (show.fit) {
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(dis)
    fit.y <- T.fit(rm, design = dis, step.method = step.method, 
                   min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction)
    betas <- fit.y$coefficients
  }
  for (i in 1:ncol(groups)) {
    group <- g[, i]
    if ((show.fit) && !is.null(betas)) {
      li <- c(2:6)
      a <- reg.coeffs(coefficients = betas, groups.vector = groups.vector, 
                      group = colnames(groups)[i])
      a <- c(a, rep(0, (7 - length(a))))
      curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) + 
              a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), from = min(time), 
            to = max(time), col = color1[i], add = TRUE, 
            lty = li[i])
    }
    if (show.lines) {
      lx <- abcissa[group != 0]
      ly <- y[group != 0]
      ord <- order(lx)
      lxo <- lx[ord]
      lyo <- ly[ord]
      lines(lxo, lyo, col = color1[i])
    }
  }
  op <- par(bg = "white")
  if (legend) 
    legend(max(abcissa, na.rm = TRUE) * 1.02, ylim[1], legend = codeg, 
           text.col = color1, col = color1, cex = cexlab, lty = 1, 
           yjust = 0)
  par(op)
}

PlotTwoGroups = function (data, data2, edesign = NULL, time = edesign[, 1], groups = edesign[, c(3:ncol(edesign))], 
                          repvect = edesign[, 2], show.fit = FALSE, x.labels = NULL,
                          dis = NULL, step.method = "backward", min.obs = 2, alfa = 0.05, 
                          nvar.correction = FALSE, summary.mode = "median", show.lines = TRUE, 
                          groups.vector = NULL, xlab = "Time", ylab = "Expression value", 
                          cex.xaxis = 1, ylim = NULL, main = NULL, cexlab = 0.8, legend = TRUE, 
                          sub = NULL, item = NULL, log2scale = TRUE) 
{
  time.change = 1:length(unique(time))
  names(time.change) = unique(time)
  time = as.numeric(time.change[as.character(time)])
  ##for data
  if (!is.vector(data)) {
    if (summary.mode == "representative") {
      distances <- apply(as.matrix(dist(data, diag = TRUE, upper = TRUE)), 1, sum)
      representative <- names(distances)[distances == min(distances)]
      yy <- as.numeric(data[rownames(data) == representative, ])
      sub <- paste("Representative:", representative)
    } else if (summary.mode == "median") {
      yy <- apply(as.matrix(data), 2, median, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Median profile of", nrow(data), item, sep = " ")
      }
    } else if (summary.mode == "mean") {
      yy <- apply(as.matrix(data), 2, mean, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Mean profile of", nrow(data), item, sep = " ")
      }
    } else stop("not valid summary.mode")
    if (dim(data)[1] == 1) {
      sub <- rownames(data)
    }
  } else if (length(data) != 0) {
    yy <- as.numeric(data)
    sub <- rownames(data)
  } else stop("empty data")
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  reps <- i.rank(repvect)
  y <- vector(mode = "numeric", length = length(unique(reps)))
  x <- vector(mode = "numeric", length = length(unique(reps)))
  g <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y)) {
    y[k] <- mean(yy[reps == k], na.rm = TRUE)
    x[k] <- mean(time[reps == k])
    for (j in 1:ncol) {
      g[k, j] <- mean(groups[reps == k, j])
    }
  }
  
  ##for data2
  if (!is.vector(data2)) {
    if (summary.mode == "representative") {
      distances <- apply(as.matrix(dist(data2, diag = TRUE, upper = TRUE)), 1, sum)
      representative <- names(distances)[distances == min(distances)]
      yy2 <- as.numeric(data2[rownames(data2) == representative, ])
      sub <- paste("Representative:", representative)
    } else if (summary.mode == "median") {
      yy2 <- apply(as.matrix(data2), 2, median, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Median profile of", nrow(data2), item, sep = " ")
      }
    } else if (summary.mode == "mean") {
      yy2 <- apply(as.matrix(data2), 2, mean, na.rm = TRUE)
      if (is.null(sub)) {
        sub <- paste("Mean profile of", nrow(data2), item, sep = " ")
      }
    } else stop("not valid summary.mode")
    if (dim(data2)[1] == 1) {
      sub <- rownames(data2)
    }
  } else if (length(data2) != 0) {
    yy2 <- as.numeric(data2)
    sub <- rownames(data2)
  } else stop("empty data2")
  if (is.null(ncol(groups))) {
    ncol = 1
    legend = FALSE
    codeg = "group"
  } else {
    ncol = ncol(groups)
    codeg <- as.character(colnames(groups))
  }
  reps <- i.rank(repvect)
  y2 <- vector(mode = "numeric", length = length(unique(reps)))
  x2 <- vector(mode = "numeric", length = length(unique(reps)))
  g2 <- matrix(nrow = length(unique(reps)), ncol = ncol)
  for (k in 1:length(y2)) {
    y2[k] <- mean(yy2[reps == k], na.rm = TRUE)
    x2[k] <- mean(time[reps == k])
    for (j in 1:ncol) {
      g2[k, j] <- mean(groups[reps == k, j])
    }
  }
  if(log2scale){
    yy = log2(yy+1)
    yy2 = log2(yy2+1)
    y = log2(y+1)
    y2 = log2(y2+1)
  }
  zz=c(yy,yy2)
  ylim = c(min(as.numeric(zz), na.rm = TRUE), max(as.numeric(zz), na.rm = TRUE)) #if (is.null(ylim)) 
  abcissa <- x
  xlim = c(min(abcissa, na.rm = TRUE), max(abcissa, na.rm = TRUE)) # * 1.3
  color1 <- as.numeric(sort(factor(colnames(groups)))) + 1
  color2 <- groups
  for (j in 1:ncol) {
    color2[, j] <- color2[, j] * j
  }
  color2 <- as.vector(apply(color2, 1, sum) + 1)
  plot(x = time, y = yy, pch = 21, xlab = xlab, ylab = ylab, 
       xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim, 
       cex = cexlab, col = "red") #
  points(x = time, y = yy2, pch = 22, xlab = xlab, ylab = ylab, 
         xaxt = "n", main = main, sub = sub, ylim = ylim, xlim = xlim, 
         cex = cexlab, col = "darkblue")
  if (is.null(x.labels)){
    axis(1, at = unique(abcissa), labels = FALSE, cex.axis = cex.xaxis)
    text(unique(abcissa),par("usr")[3]-0.2, adj = 1, labels = unique(abcissa), srt=45, xpd = NA, cex = cexlab)
  } else {
    axis(1, at = unique(abcissa), labels = FALSE, cex.axis = cex.xaxis)
    text(unique(abcissa),par("usr")[3]-0.2, adj = 1, labels = x.labels, srt=45, xpd = NA, cex = cexlab)
  }
  if (show.fit) {
    rm <- matrix(yy, nrow = 1, ncol = length(yy))
    rownames(rm) <- c("ratio medio")
    colnames(rm) <- rownames(dis)
    fit.y <- T.fit(rm, design = dis, step.method = step.method, 
                   min.obs = min.obs, alfa = alfa, nvar.correction = nvar.correction)
    betas <- fit.y$coefficients
  }
  for (i in 1:ncol(groups)) {
    group <- g[, i]
    if ((show.fit) && !is.null(betas)) {
      li <- c(2:6)
      a <- reg.coeffs(coefficients = betas, groups.vector = groups.vector, 
                      group = colnames(groups)[i])
      a <- c(a, rep(0, (7 - length(a))))
      curve(a[1] + a[2] * x + a[3] * (x^2) + a[4] * (x^3) + 
              a[5] * (x^4) + a[6] * (x^5) + a[7] * (x^5), from = min(time), 
            to = max(time), col = color1[i], add = TRUE, 
            lty = li[i])
    }
    if (show.lines) {
      lx <- abcissa[group != 0]
      ly <- y[group != 0]
      ord <- order(lx)
      lxo <- lx[ord]
      lyo <- ly[ord]
      lines(lxo, lyo, col = "red")
      lx2 <- abcissa[group != 0]
      ly2 <- y2[group != 0]
      ord2 <- order(lx2)
      lxo2 <- lx2[ord2]
      lyo2 <- ly2[ord2]
      lines(lxo2, lyo2, col = "darkblue")
    }
  }
  op <- par(bg = "white")
  if (legend) 
    legend(max(abcissa, na.rm = TRUE) * 0.8, ylim[1], legend = c("Major","Minor"), 
           text.col = c("red","darkblue"), col = c("red","darkblue"), cex = cexlab, lty = 1, 
           yjust = 0) # * 1.02
  par(op)
}

get.siggenes = function (tstep, rsq = 0.7, add.IDs = FALSE, IDs = NULL, matchID.col = 1, 
                         only.names = FALSE, vars = c("all", "each", "groups"), significant.intercept = "dummy", 
                         groups.vector = NULL, trat.repl.spots = "none", index = IDs[, (matchID.col + 1)], 
                         match = IDs[, matchID.col], r = 0.7) 
{
  dis <- tstep$dis
  edesign <- tstep$edesign
  groups.vector <- tstep$groups.vector
  not.all.empty <- function(x) (is.element(FALSE, x == " "))
  indep <- strsplit(colnames(dis)[grep("x", colnames(dis))[1]], "x")[[1]][1]
  if (any(tstep$sol[, 2] > rsq)) {
    sig.pvalues <- tstep$sol[which(tstep$sol[, 2] > rsq),]
    sig.profiles <- tstep$sig.profiles[which(tstep$sol[,2] > rsq), ]
    coefficients <- tstep$coefficients[which(tstep$sol[,2] > rsq), ]
    group.coeffs <- tstep$group.coeffs[which(tstep$sol[,2] > rsq), ]
    if (vars == "all") {
      sigs <- sig.profiles
      summary <- rownames(sig.profiles)
      if (only.names) {
        sigs <- rownames(sig.profiles)
      }
      if (add.IDs) {
        row.names <- rownames(sig.profiles)
        ids <- IDs[match(rownames(sig.profiles), IDs[,matchID.col]), ]
        if (!only.names) {
          sigs <- cbind(ids, sigs)
        }else {
          sigs <- ids
        }
        rownames(sigs) <- row.names
      }
      coeffs <- coefficients
      gc <- group.coeffs
      ps <- sig.pvalues
      sig.genes <- list(sigs, coeffs, gc, ps, nrow(sig.profiles), 
                        edesign, groups.vector)
      names(sig.genes) <- c("sig.profiles", "coefficients", 
                            "group.coeffs", "sig.pvalues", "g", "edesign", 
                            "groups.vector")
    }else if (vars == "each") {
      sig.genes <- as.list(paste("var", c("independ", colnames(dis)), sep = "."))
      summary <- matrix(" ", ncol = ncol(dis) + 1, nrow = nrow(sig.profiles))
      colnames(summary) <- c("independ", colnames(dis))
      for (i in 1:ncol(summary)) {
        sigs <- sig.profiles[which(!is.na(sig.pvalues[, (2 + i)])), ]
        coeffs <- coefficients[which(!is.na(sig.pvalues[, (2 + i)])), ]
        gc <- group.coeffs[which(!is.na(sig.pvalues[, (2 + i)])), ]
        ps <- sig.pvalues[which(!is.na(sig.pvalues[, (2 + i)])), ]
        if (nrow(sigs) > 0) names.sigs <- rownames(sigs) else names.sigs <- NULL
        summary[, i] <- c(names.sigs, rep(" ", nrow(sig.profiles) - nrow(sigs)))
        sig.genes[[i]] <- list(sigs, coeffs, gc, ps, nrow(sigs), edesign, groups.vector)
        names(sig.genes[[i]]) <- c("sig.profiles", "coefficients", 
                                   "group.coeffs", "sig.pvalues", "g", "edesign", 
                                   "groups.vector")
      }
      names(sig.genes) <- c("independ", colnames(dis))
      summary <- as.data.frame(summary[apply(summary, 1, 
                                             not.all.empty), ])
    } else if (vars == "groups") {
      if (is.null(groups.vector)) {
        if (is.null(tstep$groups.vector)) {
          stop("groups.vector is missing")
        }
        else {
          groups.vector <- tstep$groups.vector
        }
      }
      group <- unique(groups.vector)
      summary <- matrix(" ", ncol = length(group), nrow = nrow(sig.profiles))
      colnames(summary) <- group
      sig.genes <- as.list(group)
      if (significant.intercept == "all") {
        selc <- c(1:length(groups.vector))
      } else if (significant.intercept == "dummy") {
        selc <- c(2:length(groups.vector))
      } else if (significant.intercept == "none") {
        selc <- grep(indep, colnames(tstep$coefficients))
      } else stop("invalid significant.intercept value, must be one of: all, dummy, none")
      for (i in 1:length(group)) {
        group.sig <- sig.pvalues[, grep("p.valor", colnames(sig.pvalues))]
        cnames1 <- colnames(group.sig)
        group.sig <- as.data.frame(group.sig[, selc])
        if (length(selc) == 1) {
          colnames(group.sig) <- cnames1[selc]
        }
        cnames2 <- colnames(group.sig)
        group.sig <- as.data.frame(group.sig[, groups.vector[selc] == group[i]])
        if (ncol(group.sig) == 1) {
          group.sig <- as.data.frame(group.sig)
          colnames(group.sig) <- cnames2[groups.vector[selc] == group[i]]
        }
        print(length(which(apply(group.sig, 1, function(x) { any(!is.na(x)) }))))
        ps <- sig.pvalues[which(apply(group.sig, 1, function(x) { any(!is.na(x)) })), ]
        sigs <- sig.profiles[which(apply(group.sig, 1, function(x) { any(!is.na(x)) })), ]
        coeffs <- coefficients[which(apply(group.sig, 1, function(x) { any(!is.na(x)) })), ]
        gc <- group.coeffs[which(apply(group.sig, 1, function(x) { any(!is.na(x)) })), ]
        if (nrow(sigs) > 0) 
          names.sigs <- rownames(sigs)
        else names.sigs <- NULL
        summary[, i] <- c(names.sigs, rep(" ", nrow(sig.profiles) - 
                                            nrow(sigs)))
        sig.genes[[i]] <- list(sigs, coeffs, gc, ps, 
                               nrow(ps), edesign, groups.vector)
        names(sig.genes[[i]]) <- c("sig.profiles", "coefficients", 
                                   "group.coeffs", "sig.pvalues", "g", "edesign", 
                                   "groups.vector")
      }
      names(sig.genes) <- unique(groups.vector)
      if (nrow(summary) > 1) 
        summary <- as.data.frame(summary[apply(summary,1, not.all.empty), ])
    } else stop("invalid vars value, must be one of: all, each, groups")
    if (trat.repl.spots == "average") {
      if (vars != "all") {
        for (i in 1:length(sig.genes)) {
          sig.genes[[i]][[1]] <- average.rows(sig.genes[[i]][[1]], 
                                              index = index, match = match, r = r)
          for (j in c(2:3)) {
            sig.genes[[i]][[j]] <- average.rows(sig.genes[[i]][[j]], 
                                                index = index, match = match, r = -1)
            sig.genes[[i]][[j]] <- sig.genes[[i]][[j]][is.element(sig.genes[[i]][[j]], 
                                                                  sig.genes[[i]][[1]]), ]
          }
        }
      } else {
        sig.genes[[1]] <- average.rows(sig.genes[[1]], 
                                       index = index, match = match, r = r)
        for (j in c(2:3)) {
          sig.genes[[j]] <- average.rows(sig.genes[[j]], 
                                         index = index, match = match, r = -1)
          sig.genes[[j]] <- sig.genes[[j]][is.element(sig.genes[[j]], 
                                                      sig.genes[[1]]), ]
        }
      }
    }
    sig.genes2 <- sig.genes
    if (only.names && vars != "all") {
      for (i in 1:length(sig.genes)) {
        if (!is.null(dim(sig.genes[[i]][[1]]))) {
          sig.genes[[i]][[1]] <- rownames(sig.genes[[i]][[1]])
        }
      }
    }
    if (add.IDs && vars != "all") {
      for (i in 1:length(sig.genes)) {
        if (nrow(sig.genes2[[i]][[1]]) > 1) {
          row.names <- rownames(sig.genes2[[i]][[1]])
          if (trat.repl.spots == "none") {
            ids <- IDs[match(rownames(sig.genes2[[i]][[1]]), 
                             IDs[, matchID.col]), ]
          } else {
            stop("function parameters no compatible (add.IDs, trat.repl.spots)")
          }
          if (!only.names) {
            sig.genes[[i]][[1]] <- cbind(ids, sig.genes[[i]][[1]])
          }
          else sig.genes[[i]][[1]] <- ids
          rownames(sig.genes[[i]][[1]]) <- row.names
        }
      }
    }
  } else {
    sig.genes <- NULL
    summary <- c("no significant genes")
    print("no significant genes")
  }
  output <- list(sig.genes, summary)
  names(output) <- c("sig.genes", "summary")
  output
}
