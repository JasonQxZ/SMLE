
<!-- PROJECT LOGO -->
<br />
<p align="right">
  <a href="https://github.com/JasonQxZ/SMLE">
    <img src="https://user-images.githubusercontent.com/46462586/117344253-2dee5c00-ae73-11eb-8628-46b7967656f7.png" alt="Logo" width="100">
  </a>

  <h1 align="left">SMLE R Package</h1>

  <p align="left">
   <img src="http://www.r-pkg.org/badges/version/SMLE">    <img src="https://cranlogs.r-pkg.org/badges/grand-total/SMLE">
  </p>

  <p align="left">
    Joint Feature Screening via Sparse MLE
    <br />
    <a href="https://github.com//JasonQxZ/SMLE"><strong>Explore the docs »</strong></a>
  

<!-- TABLE OF CONTENTS -->

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>


## Getting Started
------------

### Prerequisites

R softeware and R package glmnet


### Installation

    # Install package from CRAN:
    install.packages("SMLE")


## Usage 

    Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
    fit<-SMLE(Data$Y, Data$X, k=9, family = "gaussian")
  
For more examples, please refer to the [Vignette](http://https://github.com/JasonQxZ/SMLE/blob/main/vignette/vignette.md)
## License

Distributed under the GPL-3 License.

## Contact
Qianxiang Zang - qzang023@uottawa.ca
[![LinkedIn][linkedin-shield]][linkedin-url]


<!-- MARKDOWN LINKS & IMAGES -->
[forks-shield]: https://img.shields.io/github/forks/JasonQxZ/SMLE.svg?style=for-the-badge
[forks-url]: https://github.com/JasonQxZ/SMLEnetwork/members
[stars-shield]: https://img.shields.io/github/stars/JasonQxZ/SMLE.svg?style=for-the-badge
[stars-url]: https://github.com/asonQxZ/SMLE/stargazers
[issues-shield]: https://img.shields.io/github/issues/asonQxZ/SMLE.svg?style=for-the-badge
[issues-url]: https://github.com/asonQxZ/SMLE/issues
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/www.linkedin.com/in/qianxiangzang

