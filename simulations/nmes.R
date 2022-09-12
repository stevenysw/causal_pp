smoke <- read.csv("nmesdata.txt", header = T)
smoke <- subset(smoke, TOTALEXP > 0)
smoke <- subset(smoke, select = c(TOTALEXP, packyears, LASTAGE, AGESMOKE,
                    MALE, RACE3, marital, educate,
                     SREGION, POVSTALB, beltuse, yearsince))
smoke <- na.omit(smoke)
smoke$weight <- nmes_data[row.names(smoke),]$HSQACCWT
smoke <- na.omit(smoke)

X = as.matrix(smoke[,-c(1,2)])
y = log(smoke$TOTALEXP)
z = as.numeric(smoke$packyears > 17)

rpart.plot(tau_cart$model)
