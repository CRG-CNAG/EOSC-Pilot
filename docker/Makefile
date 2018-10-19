gatk1.2:
	docker build -t cbcrg/eosc-pilot:1.2 . -f Dockerfile-1.2

gonl3.4:
	docker build -t cbcrg/eosc-pilot:gonl3.4 . -f Dockerfile.gonl3.4

gatk4:
	docker build -t cbcrg/eosc-pilot:gatk4 . -f Dockerfile.gatk4

travis:
	docker build -f Dockerfile.travis -t 885800555707.dkr.ecr.eu-west-1.amazonaws.com/cbcrg/eosc_with_gatk3.7 .
	docker push 885800555707.dkr.ecr.eu-west-1.amazonaws.com/cbcrg/eosc_with_gatk3.7