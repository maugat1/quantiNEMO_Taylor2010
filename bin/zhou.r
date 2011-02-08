#!/usr/bin/R -q --no-save

library("Hmisc");

################################################################################
# settings
################################################################################
files = c(
        'results/zhou_prop_N25_s1e-1/patch_inbr.csv',
        'results/zhou_prop_N100_s1e-1/patch_inbr.csv',
        'results/zhou_prop_N25_s1e-2/patch_inbr.csv',
        'results/zhou_prop_N100_s1e-2/patch_inbr.csv'
        );
models = c(
        'propagule',
        'propagule',
        'propagule',
        'propagule'
        );
popsizes = c(
        25,
        100,
        25,
        100
        );
selections = c(
        0.1,
        0.1,
        0.01,
        0.01
        );

parameters = paste(models, paste("N=", popsizes, sep=''),
           paste("s=", selections, sep=''), sep=", ");
#colors = rainbow(length(files)+7);
colors = c(
        "#FF0000",
        "#00AA55",
        "#0000FF",
        "#AA00AA"
        );
#linetype = c(1:length(files));
linetype = c(1, 1, 1, 1);
plotchar = seq(18, 18 + length(files), 1);

outfiles = c(
        "prop_inbr_depr.pdf",
        "prop_gen_load.pdf"
        );
width = 6;
height = 6;

scale_sdev = 10;

################################################################################
# get the data ranges
################################################################################
ages.range = c(0, 0);
inbr.range = c(0, 0);
load.range = c(0, 0);
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.inbr = patchtab$inbrDepr;
    patch.load = 1 - patchtab$meanW;
    # get the mean and variance for each coln age
    next_ind = 1
    ages = c();
    inbr.mean = c();
    inbr.sdev = c();
    load.mean = c();
    load.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.inbr);
        if(any(is_j)) {
            ages[next_ind] = j;
            inbr.mean[next_ind] = mean(patch.inbr[is_j]);
            inbr.sdev[next_ind] = sqrt(mean((patch.inbr[is_j] -
                            inbr.mean[next_ind])**2)) / scale_sdev;
            load.mean[next_ind] = mean(patch.load[is_j]);
            load.sdev[next_ind] = sqrt(mean((patch.load[is_j] -
                            load.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    ages.range[1] = min(ages.range[1], min(ages));
    ages.range[2] = max(ages.range[2], max(ages));
    inbr.range[1] = min(inbr.range[1], min(inbr.mean - inbr.sdev));
    inbr.range[2] = max(inbr.range[2], max(inbr.mean + inbr.sdev));
    load.range[1] = min(load.range[1], min(load.mean - load.sdev));
    load.range[2] = max(load.range[2], max(load.mean + load.sdev));
}

################################################################################
# plot the inbreeding depression
################################################################################
pdf(outfiles[1], width=width, height=height);
plot(ages.range, inbr.range, type='n', xlab="Patch age from colonization",
        ylab="Inbreeding depression");
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.inbr = patchtab$inbrDepr;
    # get the mean and variance for each coln age
    next_ind = 1;
    ages = c();
    inbr.mean = c();
    inbr.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.inbr);
        if(any(is_j)) {
            ages[next_ind] = j;
            inbr.mean[next_ind] = mean(patch.inbr[is_j]);
            inbr.sdev[next_ind] = sqrt(mean((patch.inbr[is_j] -
                            inbr.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    lines(ages, inbr.mean, type='b', lwd=1.5, lty=linetype[i], pch=plotchar[i],
            col=colors[i]);
    errbar(ages, inbr.mean, inbr.mean + inbr.sdev, inbr.mean - inbr.sdev,
            add=TRUE, col=colors[i]);
}
legend(ages.range[1], inbr.range[2], parameters, col=colors,
        pch=plotchar, lty=linetype);
dev.off();

################################################################################
# plot the genetic load
################################################################################
pdf(outfiles[2], width=width, height=height);
plot(ages.range, load.range, type='n', xlab="Patch age from colonization",
        ylab="Genetic load");
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.load = 1 - patchtab$meanW;
    # get the mean and variance for each coln age
    next_ind = 1;
    ages = c();
    load.mean = c();
    load.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.load);
        if(any(is_j)) {
            ages[next_ind] = j;
            load.mean[next_ind] = mean(patch.load[is_j]);
            load.sdev[next_ind] = sqrt(mean((patch.load[is_j] -
                            load.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    lines(ages, load.mean, type='b', lwd=1.5, lty=linetype[i], pch=plotchar[i],
            col=colors[i]);
    errbar(ages, load.mean, load.mean + load.sdev, load.mean - load.sdev,
            add=TRUE, col=colors[i]);
}
legend(ages.range[1], load.range[1] + 0.03, parameters, col=colors,
        pch=plotchar, lty=linetype);
dev.off();



################################################################################
# settings
################################################################################
files = c(
        'results/zhou_migr_N25_s1e-1/patch_inbr.csv',
        'results/zhou_migr_N100_s1e-1/patch_inbr.csv',
        'results/zhou_migr_N25_s1e-2/patch_inbr.csv',
        'results/zhou_migr_N100_s1e-2/patch_inbr.csv'
        );
models = c(
        'migrant',
        'migrant',
        'migrant',
        'migrant'
        );
popsizes = c(
        25,
        100,
        25,
        100
        );
selections = c(
        0.1,
        0.1,
        0.01,
        0.01
        );

parameters = paste(models, paste("N=", popsizes, sep=''),
           paste("s=", selections, sep=''), sep=", ");
#colors = rainbow(length(files)+7);
colors = c(
        "#FF0000",
        "#00AA55",
        "#0000FF",
        "#AA00AA"
        );
#linetype = c(1:length(files));
linetype = c(1, 1, 1, 1);
plotchar = seq(18, 18 + length(files), 1);

outfiles = c(
        "migr_inbr_depr.pdf",
        "migr_gen_load.pdf"
        );

################################################################################
# get the data ranges
################################################################################
ages.range = c(0, 0);
inbr.range = c(0, 0);
load.range = c(0, 0);
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.inbr = patchtab$inbrDepr;
    patch.load = 1 - patchtab$meanW;
    # get the mean and variance for each coln age
    next_ind = 1
    ages = c();
    inbr.mean = c();
    inbr.sdev = c();
    load.mean = c();
    load.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.inbr);
        if(any(is_j)) {
            ages[next_ind] = j;
            inbr.mean[next_ind] = mean(patch.inbr[is_j]);
            inbr.sdev[next_ind] = sqrt(mean((patch.inbr[is_j] -
                            inbr.mean[next_ind])**2)) / scale_sdev;
            load.mean[next_ind] = mean(patch.load[is_j]);
            load.sdev[next_ind] = sqrt(mean((patch.load[is_j] -
                            load.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    ages.range[1] = min(ages.range[1], min(ages));
    ages.range[2] = max(ages.range[2], max(ages));
    inbr.range[1] = min(inbr.range[1], min(inbr.mean - inbr.sdev));
    inbr.range[2] = max(inbr.range[2], max(inbr.mean + inbr.sdev));
    load.range[1] = min(load.range[1], min(load.mean - load.sdev));
    load.range[2] = max(load.range[2], max(load.mean + load.sdev));
}

################################################################################
# plot the inbreeding depression
################################################################################
pdf(outfiles[1], width=width, height=height);
plot(ages.range, inbr.range, type='n', xlab="Patch age from colonization",
        ylab="Inbreeding depression");
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.inbr = patchtab$inbrDepr;
    # get the mean and variance for each coln age
    next_ind = 1;
    ages = c();
    inbr.mean = c();
    inbr.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.inbr);
        if(any(is_j)) {
            ages[next_ind] = j;
            inbr.mean[next_ind] = mean(patch.inbr[is_j]);
            inbr.sdev[next_ind] = sqrt(mean((patch.inbr[is_j] -
                            inbr.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    lines(ages, inbr.mean, type='b', lwd=1.5, lty=linetype[i], pch=plotchar[i],
            col=colors[i]);
    errbar(ages, inbr.mean, inbr.mean + inbr.sdev, inbr.mean - inbr.sdev,
            add=TRUE, col=colors[i]);
}
legend(ages.range[1], inbr.range[1] + 0.03, parameters, col=colors,
        pch=plotchar, lty=linetype);
dev.off();

################################################################################
# plot the genetic load
################################################################################
pdf(outfiles[2], width=width, height=height);
plot(ages.range, load.range, type='n', xlab="Patch age from colonization",
        ylab="Genetic load");
for(i in 1:length(files)) {
    infile = files[i];
    # read the data
    patchtab = read.table(infile, header=TRUE, sep=",", na.strings="NaN",
            strip.white=TRUE, blank.lines.skip=TRUE);
    patch.age = patchtab$colnAge;
    patch.load = 1 - patchtab$meanW;
    # get the mean and variance for each coln age
    next_ind = 1;
    ages = c();
    load.mean = c();
    load.sdev = c();
    last_j = -1;
    for(j in (1:max(patch.age))[c(FALSE, FALSE, FALSE, FALSE, TRUE)]) {
        is_j = (patch.age > last_j) & (patch.age <= j) & !is.na(patch.load);
        if(any(is_j)) {
            ages[next_ind] = j;
            load.mean[next_ind] = mean(patch.load[is_j]);
            load.sdev[next_ind] = sqrt(mean((patch.load[is_j] -
                            load.mean[next_ind])**2)) / scale_sdev;
            next_ind = next_ind + 1;
        }
        last_j = j;
    }
    lines(ages, load.mean, type='b', lwd=1.5, lty=linetype[i], pch=plotchar[i],
            col=colors[i]);
    errbar(ages, load.mean, load.mean + load.sdev, load.mean - load.sdev,
            add=TRUE, col=colors[i]);
}
legend(ages.range[1], load.range[1] + 0.03, parameters, col=colors,
        pch=plotchar, lty=linetype);
dev.off();
