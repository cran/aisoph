\name{aisoph}
\alias{aisoph}
\title{Fit Additive Isotonic Proportional Hazards Model}
\description{Nonparametric estimation of additive isotonic covariate effects for proportional hazards model.}
\usage{
  aisoph(time, status, z1, z2, w, shape, K1, K2, maxiter, eps)
}
\arguments{
  \item{time}{survival time. It must be greater than 0.}
  \item{status}{censoring indication. It must be 0 or 1.}
  \item{z1}{First covariate under order-restriction.}
  \item{z2}{Second covariate under-order restriction.}
  \item{w}{Additional covariates. This argument is optional}
  \item{shape}{Directions for \eqn{ z1 } and \eqn{ z2 }, "increasing" or "decreasing".}
  \item{K1}{anchor constraint for \eqn{ z1 } (default is 0).}
  \item{K2}{anchor constraint for \eqn{ z2 } (default is 0).}
  \item{maxiter}{maximum number of iteration (default is 10^5).}
  \item{eps}{stopping convergence criteria (default is 10^-3).}
}
\details{The aisoph function allows to analyze additive isotonic proportional hazards model, which is defined as
\deqn{ \lambda(t|z1, z2, w)=\lambda0(t)exp(\psi1(z1)+\psi2(z2)+\beta w), }
where \eqn{ \lambda0 } is a unspecified baseline hazard function, \eqn{ \psi1 } and \eqn{ \psi2 } are monotone increasing (or decreasing) functions in \eqn{z1} and \eqn{z2}, respectively, \eqn{w} is additional covariate(s), and \eqn{ \beta } is a regression paramter(s). If \eqn{ w } is omitted in the formulation above, \eqn{ \psi1 } and \eqn{ \psi2 } are estimated as right continuous increasing (or left continuous decreasing) step functions. Otherwise, \eqn{ \psi1, \psi2 } and \eqn{ w } are estimated.

For the anchor constraint, one point has to be fixed with \eqn{\psi1(K1)=0} to solve the identifiability problem, e.g. \eqn{\lambda0(t)exp(\psi1(z1)+\psi2(z2)+\beta w)=(\lambda0(t)exp(-c))exp(\psi1(z)+c+\psi2(z2)+\beta w)} for any constant \eqn{c}. Similarly, \eqn{\psi2(K2)=0}. \eqn{K1} and \eqn{K2} are called anchor points. By default, we set \eqn{K1} and \eqn{K2} as medians of \eqn{z1} s and \eqn{z2} values, respectively. The choise of anchor points are less importants because hazars ratios are not affected by anchor points.}
\value{A list of class isoph:
\item{iso1}{data.frame with estimated \eqn{z1}, \eqn{\psi1}, and \eqn{\exp(\psi1)}.}
\item{iso2}{data.frame with estimated \eqn{z2}, \eqn{\psi2}, and \eqn{\exp(\psi2)}.}
\item{est}{estimated \eqn{\beta}, and \eqn{\exp(\beta)}}
\item{conv}{Algorithm convergence status.}
\item{shape1}{shape-constrain for \eqn{\psi1}.}
\item{shape2}{shape-constrain for \eqn{\psi2}.}
}

\references{Yunro Chung, Anastasia Ivanova, Jason P. Fine, Additive isotonic proportional hazards models (working in progress).}
\author{Yunro Chung [aut, cre]}

\examples{
#require(survival)
#require(Iso)

###
# 1. time-independent covariate with monotone increasing effect
###
# 1.1. create a test data set 1
time=  c(1, 6, 3, 6, 7, 8, 1, 4, 0, 2, 1, 5, 8, 7, 4)
status=c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
z1=     c(3, 1, 2, 4, 8, 3, 3, 4, 1, 9, 4, 2, 2, 8, 5)
z2=     c(1, 3, 5, 6, 1, 7, 6, 8, 3, 4, 8, 8, 5, 2, 3)

# 1.2. Fit isotonic proportional hazards model
res1 = aisoph(time=time, status=status, z1=z1, z2=z2, shape="increasing")

# 1.3. print result
res1

#1.4. plot
plot(res1)

###
# 2. time-independent covariate with monotone increasing effect
###
# 2.1. create a test data set 1
time=  c(0,4,8,9,5,6,9,8,2,7,4,2,6,2,5,9,4,3,8,2)
status=c(0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1)
z1=    c(3,2,1,1,3,1,8,4,3,6,2,9,9,0,7,7,2,3,4,6)
z2=    c(3,6,9,9,4,3,9,8,4,7,2,3,1,3,7,0,1,6,4,1)
trt=   c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)

# 2.2. Fit isotonic proportional hazards model
res2 = aisoph(time=time, status=status, z1=z1, z2=z2, w=trt, shape="increasing")

# 2.3. print result
res2

#2.4. plot
plot(res2)
}

\concept{Isotonic regression}
\concept{Survival analysis}
\concept{Nonparametric regression}
\concept{Order-restricted inference}
