# saaga-api

The Species-agnostic Automated Gas Analyzer (SAAGA) project aims to automate the detection and characterization of chemical compounds in a complex chemical mixture of gas phase through experimental rotational spectroscopy and computational tools. A database of spectroscopic data serves as the foundation of the automation pipeline for assigning spectral lines to species. While there are existing databases available for use, we developed our custom database, named SAAGAdb, and an application programming interface (API) to access the database to fulfil the needs of SAAGA.

SAAGAdb is designed to store structured, high quality spectroscopic data of all species not limited to astrochemically relevant ones, enabling convenient data manipulation, integration into future automation pipeline, deployment, and maintenance. We implemented software development best practices, including software development life cycle, continuous integration/continuous delivery, and version control, to develop a PostgreSQL database with a Python API built on Django with RDKit integration. SAAGAdb and its API are up and running on a physical server. The RESTful API allows CRUD operations of the database, in which admins can perform the full CRUD operations while other users can read and query the database. The product passed all unit tests and was successfully seeded with data. With the flexibility provided by the Django framework as well as detailed documentation of the software, SAAGAdb and its API can be easily improved and expanded in the future to suit the needs of the SAAGA project.