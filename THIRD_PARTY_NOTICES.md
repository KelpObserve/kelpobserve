# Third-party notices

The code and documents of this repository are licensed under the Apache License 2.0 (`LICENSE`).

## Data

| Source | Use | Attribution and terms |
|---|---|---|
| Copernicus Sentinel-2, Level-1C and Level-2A | Imagery of the 2026 campaigns and of the winter scenes; Map 4 of the technical report | Contains modified Copernicus Sentinel data 2025–2026. |
| Copernicus Marine Service, IBI sea level: dataset `cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m` | Water level at each acquisition | Generated using E.U. Copernicus Marine Service Information. |
| Copernicus Marine Service, Atlantic ocean colour: product `OCEANCOLOUR_ATL_BGC_L3_MY_009_113` ([doi:10.48670/moi-00286](https://doi.org/10.48670/moi-00286)) | KD490 for the exploratory attenuation module | Generated using E.U. Copernicus Marine Service Information. |
| EMODnet Bathymetry: the high-resolution lidar model `590_HR_Lidar_Sul` and the EMODnet DTM | Depth. `benchmark/inputs/emodnet/` holds both, made on 7 September 2026: cut to the Sesimbra area with land masked out, and for the DTM, read through EMODnet's coverage service, depth made positive down | EMODnet Bathymetry Consortium. The DTM was created by EMODnet (https://emodnet.ec.europa.eu/en/), is owned by the EU and is licensed under CC BY 4.0. `590_HR_Lidar_Sul`: lidar survey of 2011 by DGT and APA, distributed by Instituto Hidrográfico (EDMO 590), with the use limitation "public" in its EMODnet metadata. Not for navigation. |
| EMODnet Seabed Habitats: survey map PT010000 (LIFE Biomares and IPMA, 2013) and EUSeaMap 2025 | Habitat context (technical report, section 2.6) | EMODnet Seabed Habitats; PT010000 from LIFE Biomares and IPMA (2013). |
| Pléiades Neo, supplied through ESA Third Party Missions project PP0110832 | Imagery of 8 July and 3 September 2026 | Pléiades Neo © Airbus DS (2026). No Pléiades Neo imagery is included. |

### Spectral libraries used by E5s

Only statistics derived from the spectra appear in this repository; the spectra themselves are not redistributed.

| Library | Record | Terms |
|---|---|---|
| Rowan, Gillian, Phinn, Stuart, Smart, Joanna N. and Roelfsema, Chris (2025). Hyperspectral Reflectance Measurements of Kelps, Macroalgae, and Co-ocurring Materials Collected in Australia and New Zealand. The University of Queensland. (Dataset) doi: 10.48610/c5f745e | UQ eSpace, [doi:10.48610/c5f745e](https://doi.org/10.48610/c5f745e) | Permitted re-use with acknowledgement |
| Vahtmäe et al. (2025): Baltic macroalgae and substrates | PANGAEA, [doi:10.1594/PANGAEA.971518](https://doi.org/10.1594/PANGAEA.971518) | CC BY 4.0 |
| Fuller, Martin and Asner (2024): algae of Hawaii | Zenodo, [doi:10.5281/zenodo.10826300](https://doi.org/10.5281/zenodo.10826300) | CC BY 4.0 |
| Roelfsema and Phinn (2013): Glovers Reef | PANGAEA, [doi:10.1594/PANGAEA.824861](https://doi.org/10.1594/PANGAEA.824861) | CC BY 3.0 |
| Ecotone AS: Ria de Vigo seafloor sands and sediments | Zenodo, [doi:10.5281/zenodo.13381713](https://doi.org/10.5281/zenodo.13381713) | CC BY 4.0 |
| Anstee et al. (2017): Wallis Lake | PANGAEA, [doi:10.1594/PANGAEA.873219](https://doi.org/10.1594/PANGAEA.873219) | CC BY 3.0 |

## Map backgrounds of the report figures

| Background | Figures | Credit |
|---|---|---|
| Esri World Imagery | Maps 3 and 5 of the technical report (`benchmark/docs/figures/map-3-transects.jpg` and `map-5-habitat.jpg`) | Imagery: Esri, Maxar, Earthstar Geographics and the GIS user community |
| OceanStream basemap | Maps 1 and 5 of the technical report | Basemap: OceanStream (GEBCO, EMODnet, Copernicus DEM, © OpenStreetMap contributors) |
| EMODnet Bathymetry | Map 2 of the technical report | Bathymetry: EMODnet Bathymetry Consortium |
