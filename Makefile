build:
	docker build --no-cache -t jags-rstudio -f Dockerfile .

run:
	docker run -d -p 8787:8787 -v /Users/lpatruno/dev/bayesian_inference:/home/rstudio/bayesian_inference jags-rstudio
