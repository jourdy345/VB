# Variational Bayes

Welcome to my repository for variational inference! The repository you've just stumbled upon is not for public use/exhibit; it has been created entirely for personal purpose. I don't mind people using my code or reading my summary notes but please not that there is no manual. This README file has also been written to inform MYSELF and not to forget. So you're seeing my notebook basically.

### GP mixed C++ code compiling & linking
```bash
$ g++ -c GPProbit.cpp
$ g++ -o GPPprobit GPPProbit.o -llapack -lblas -larmadillo
```
Don't forget to link the added libraries.