# 🇬🇧 Machine Learning in Drug Design (MLDD) 2021/2022

This repository contains course materials for the course "Machine Learning in Drug Design." In the `labs` directory, there are materials covering the following topics:

1. Python and machine learning basics revision. SMILES representation and molecular fingerprints. Introduction to the RDKit package.
2. Use of publicly available molecular data sets. Effective search in open databases. Creation and curation of custom data sets.
3. Protein chemistry. Popular data formats. Molecular docking.
4. Molecular graph theory. Graph neural networks. Interpretability.
5. Review of generative models in drug discovery.
6. Finding drug targets. Deep learning for proteins.

## About us

[GMUM](https://gmum.net/) (Machine Learning Research Group) is a group at the Jagiellonian University working on various aspects of machine learning, and in particular deep learning - in both fundamental and applied settings. The group is led by prof. Jacek Tabor.

Some of the research directions our group pursues include:
- generative models: efficient training and sampling; inpainting; super-resolution,
- theoretical understanding of deep learning and optimization,
- natural language processing,
- drug design and cheminformatics,
- unsupervised learning and clustering,
- computer vision and medical image analysis.

We are hosting machine learning seminars that are open to the public. You can check the schedule on [our website](https://gmum.net/seminars.html) and join online (links posted on [our Facebook](http://facebook.com/gmum.net)).
You can also add seminar info to your [Google calendar](https://calendar.google.com/calendar/u/0?cid=ZDJjcTFudnU0Y2UxNXNnODltdDc4Y3BtcTBAZ3JvdXAuY2FsZW5kYXIuZ29vZ2xlLmNvbQ).

## Environment Setup

Python will be used throughout the course. The environment setup steps are shown below:

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html) following the instructions for your operating system.
2. Download this repository: `git clone https://github.com/gmum/umwpl2021.git`.
3. Install environment from the YAML file: `conda env create -f environment.yml`

In the `environment-versions.yml` file, the exact versions of each package are listed. They may not be compatible with all operating systems.

_Important! If you would like to use your GPU to train neural networks, add `cudatoolkit` in a correct version (e.g. `cudatoolkit=10.2`) to the `environment.yml` file._

_If you plan to use our computational cluster, ask us about a configured environment on the server. We have a singularity image with all the packages installed._

## Literature

1. Rogers, D., & Hahn, M. (2010) [Extended-Connectivity Fingerprints](https://pubs.acs.org/doi/10.1021/ci100050t). *Journal of chemical information and modeling*.
2. Durant, J. L., Leland, B. A., Henry, D. R., & Nourse, J. G. (2002). [Reoptimization of MDL Keys for Use in Drug Discovery](https://pubs.acs.org/doi/10.1021/ci010132r). *Journal of chemical information and computer sciences*.
3. Deng, Z., Chuaqui, C., & Singh, J. (2004). [Structural Interaction Fingerprint (SIFt):  A Novel Method for Analyzing Three-Dimensional Protein−Ligand Binding Interactions](https://pubs.acs.org/doi/10.1021/jm030331x). *Journal of medicinal chemistry*.
4. McNutt, A. T., Francoeur, P., Aggarwal, R., Masuda, T., Meli, R., Ragoza, M., ... & Koes, D. R. (2021). [GNINA 1.0: molecular docking with deep learning](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2). *Journal of cheminformatics*.
5. Koes, D. R., Baumgartner, M. P., & Camacho, C. J. (2013). [Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise](https://pubs.acs.org/doi/abs/10.1021/ci300604z). *Journal of chemical information and modeling*.
6. Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). [Highly accurate protein structure prediction with AlphaFold](https://www.nature.com/articles/s41586-021-03819-2). *Nature*.
7. Kipf, T. N., & Welling, M. (2016). [Semi-supervised classification with graph convolutional networks](https://arxiv.org/pdf/1609.02907.pdf?fbclid=IwAR0BgJeoKHIAvPuSE9fJ0_IQOEu5l75yxyNo7PUC08RTOFlm_IIo5YmcnQM). *arXiv preprint arXiv:1609.02907*.
8. Xu, K., Hu, W., Leskovec, J., & Jegelka, S. (2018). [How powerful are graph neural networks?](https://arxiv.org/pdf/1810.00826.pdf). *arXiv preprint arXiv:1810.00826*.
9. Veličković, P., Cucurull, G., Casanova, A., Romero, A., Lio, P., & Bengio, Y. (2017). [Graph attention networks](https://arxiv.org/pdf/1710.10903.pdf). *arXiv preprint arXiv:1710.10903*.
10. Hamilton, W. L., Ying, R., & Leskovec, J. (2017, December). [Inductive representation learning on large graphs](https://proceedings.neurips.cc/paper/2017/file/5dd9db5e033da9c6fb5ba83c7a7ebea9-Paper.pdf). In *Proceedings of the 31st International Conference on Neural Information Processing Systems* (pp. 1025-1035).
11. Pope, P. E., Kolouri, S., Rostami, M., Martin, C. E., & Hoffmann, H. (2019). [Explainability methods for graph convolutional neural networks](https://openaccess.thecvf.com/content_CVPR_2019/papers/Pope_Explainability_Methods_for_Graph_Convolutional_Neural_Networks_CVPR_2019_paper.pdf). In *Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition* (pp. 10772-10781).
12. Popova, M., Isayev, O., & Tropsha, A. (2018). [Deep reinforcement learning for de novo drug design](https://www.science.org/doi/pdf/10.1126/sciadv.aap7885). *Science advances*.
13. Jin, W., Barzilay, R., & Jaakkola, T. (2018, July). [Junction tree variational autoencoder for molecular graph generation](http://proceedings.mlr.press/v80/jin18a/jin18a.pdf). In *International conference on machine learning* (pp. 2323-2332). PMLR.
14. Gómez-Bombarelli, R., Wei, J. N., Duvenaud, D., Hernández-Lobato, J. M., Sánchez-Lengeling, B., Sheberla, D., ... & Aspuru-Guzik, A. (2018). [Automatic chemical design using a data-driven continuous representation of molecules](https://pubs.acs.org/doi/pdf/10.1021/acscentsci.7b00572). *ACS central science, 4*(2), 268-276.
15. Maziarka, Ł., Danel, T., Mucha, S., Rataj, K., Tabor, J., & Jastrzębski, S. (2020). [Molecule attention transformer](https://arxiv.org/pdf/2002.08264.pdf). *arXiv preprint arXiv:2002.08264*.
16. Jiménez, J., Doerr, S., Martínez-Rosell, G., Rose, A. S., & De Fabritiis, G. (2017). [DeepSite: protein-binding site predictor using 3D-convolutional neural networks](https://academic.oup.com/bioinformatics/article/33/19/3036/3859178). *Bioinformatics, 33*(19), 3036-3042.

# 🇵🇱 Uczenie maszynowe w projektowaniu leków (UMwPL) 2021/2022

Repozytorium zawiera materiały z kursu "Uczenie maszynowe w projektowaniu leków". W folderze `labs` znajdują się materiały z ćwiczeń, które pokrywają następujące tematy:

1. Powtórka Pythona i podstaw uczenia maszynowego. Reprezentacja SMILES i fingerprinty związków chemicznych. Wprowadzenie do pakietu RDKit.
2. Korzystanie z publicznie dostępnych zbiorów danych molekuł. Efektywne przeszukiwanie otwartych baz danych. Tworzenie własnych zbiorów danych.
3. Chemia białek. Popularne formaty danych. Dokowanie molekularne.
4. Teoria grafów molekularnych. Grafowe sieci neuronowe, Interpretowalność.
5. Przegląd modeli generatywnych w odkrywaniu leków.
6. Znajdowanie nowych celów biologicznych. Uczenie głębokie dla białek.

## Instalacja Środowiska

Zajęcia będą prowadzone z użyciem języka Pythonie. Instrukcja instalacji środowiska zawierającego potrzebne paczki:

1. Instalacja [minicondy](https://docs.conda.io/en/latest/miniconda.html) zgodnie z instrukcją dla wybranego systemu.
2. Pobranie tego repozytorium: `git clone https://github.com/gmum/umwpl2021.git`.
3. Instalacja środowiska z pliku: `conda env create -f environment.yml`

W pliku `environment-versions.yml` znajdują się dokładne wersje poszczególnych paczek, ale lista może nie być kompatybilna z wszystkimi systemami.

_Ważne! Jeżeli chciałbyś użyć karty graficznej do treningu sieci neuronowych, dodaj `cudatoolkit` z odpowiednią wersją (np. `cudatoolkit=10.2`) do pliku `environment.yml`._

_Jeśli planujesz używać naszego klastra obliczeniowego, spytaj nas o skonfigurowane środowisko. Na serwerze mamy obraz singularity z zainstalownymi wszystkimi potrzebnymi paczkami._

## Projekt Zaliczeniowy

Ocena z ćwiczeń wystawiana jest na podstawie projektu. Warunkiem koniecznym uzyskania oceny pozytywnej jest uczestnictwo na zajęciach (dopuszczalna 1 nieobecność).

Projekty można będzie wykonywać w grupach 1-3 osób. W razie potrzeby możliwe jest otrzymanie dostępu do klastra obliczeniowego z GPU po złożeniu prośby do administratora. Projekt zaliczeniowy powinien posiadać:

- Publiczne repozytorium Github (lub Gitlab/Bitbucket) z licencją umożliwiającą dzielenie się kodem i wynikami projektu na stronie kursu oraz wykorzystanie kodu w przyszłych edycjach kursu (np. MIT, GNU GPL, Apache 2.0 lub własna).
- Dokumentację projektu w postaci pliku README, który zawiera między innymi: krótki opis projektu, instrukcję uruchomienia kodu wraz z listą zależności, podsumowanie wyników. Możliwe są również odniesienia do wygenerowanych plików PDF i iPython notebooków z uzupełnieniem dokumentacji, pod warunkiem że znajdują się również w repozytorium.
- Wszystkie dane lub odniesienia do źródeł danych.
- Kod umożliwiający odtworzenie kluczowych wyników projektu.

*Uwaga:* projekty, które będą kontynuowane po skończeniu semestru, np. w celu napisania publikacji, mogą mieć prywatne repozytorium udostępnione prowadzącym na tych samych zasadach.

### Ocenianie

6 i 13 XII odbędą się prezentacje projektów. Oceniane będzie zapoznanie się z tematem, zrozumienie danych i plany na pozostałą część semestru. Prezentacja będzie stanowić 15% oceny końcowej. Dla ułatwienia przygotowania prezentacji, poniżej znajdują się szczegółowe kryteria. Podane są też pytania, które mogą pomóc w przygotowaniu.

- Zrozumienie tematu (5%)
  - Co jest celem projektu?
  - Jakie jest znaczenie biologiczne/chemiczne projektu?
  - Co jest potencjalną trudnością w wykonaniu projektu?
  - Jaki jest spodziewany efekt projektu?
- Zrozumienie danych (5%)
  - jakie dane wejściowe będą użyte w projekcie?
  - skąd można pozyskać dane do projektu?
  - jaka jest struktura danych użytych w projekcie (np. fingerprinty albo grafy molekularne)?
  - wstępnie jakie problemy w dostępnych danych są widoczne?
  - jakie dodatkowe informacje (metadane) są dostępne?
- Planowana implementacja (5%)
  - Czy ten temat był już poruszany w literaturze? Jeśli tak, to jakie narzędzia są dostępne?
  - Jakie metody uczenia maszynowego będą wykorzystane w projekcie?
  - Jak zdefiniowane będzie wejście i wyjście modelu?
  - Jakie miary będą zastosowane do zmierzenia skuteczności modelu?
  - Jaki stos technologiczny będzie użyty do wykonania projektu?

Prezentacja odbywa się stacjonarnie. Jest dostęp do rzutnika, więc można przygotować pomocnicze slajdy. Prezentacja powinna zająć nie więcej niż 15 minut (najlepiej 10 minut prezentacji i 5 minut dyskusji).
