# -*- coding: utf-8 -*-
"""
CanopyMapper
============

Copyright (C) 2026  Filippo

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Script per la Processing Toolbox di QGIS che automatizza il workflow di
individuazione delle chiome (canopy) di alberi/arbusti a partire da
un'ortofoto, con calcolo automatico di area, perimetro, compattezza e,
se disponibile un DSM + DTM, dell'altezza (CHM = DSM - DTM).

WORKFLOW
--------
1. Calcolo indice di vegetazione (NDVI se è disponibile la banda NIR,
   altrimenti ExG - Excess Green - da sola ortofoto RGB)
2. Sogliatura (soglia manuale oppure automatica col metodo di Otsu)
3. (opzionale) Calcolo CHM = DSM - DTM e riallineamento sulla griglia
   dell'ortofoto; esclusione di prati/campi bassi PRIMA della
   poligonizzazione, in base a un'altezza minima
4. (opzionale) Separazione delle chiome adiacenti (watershed): trova le
   cime nel CHM (o nella forma della macchia, se niente DSM/DTM) e divide
   i gruppi fusi in poligoni singoli, uno per cima, invece di un unico
   blob da accettare o scartare in blocco
5. Poligonizzazione della maschera (o del raster etichettato dalla
   watershed) + pulizia geometrie
6. Filtro per area minima e compattezza (elimina rumore residuo)
7. (opzionale) Statistiche zonali di altezza (media, massima) per ogni
   chioma dal CHM
8. ID progressivo, arrotondamento a 2 decimali, classe dimensionale
   automatica incrociata con l'altezza (arbusto / albero piccolo / medio
   / grande / sospetto_area_ampia_bassa)

INSTALLAZIONE
--------------
Vedi il README allegato. In breve:
QGIS > Impostazioni > Opzioni utente > Sistema, oppure più semplice:
Processing Toolbox > icona Python (Script) > "Add Script to Toolbox..."
e seleziona questo file. Comparirà in Processing Toolbox sotto
"Scripts > CanopyMapper > Individua canopy da ortofoto".

Compatibile con QGIS 3.x (usa solo API core, nessuna dipendenza esterna
oltre numpy, già incluso nell'installazione standard di QGIS).
"""

import numpy as np
import heapq

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterBand,
    QgsProcessingParameterEnum,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterVectorDestination,
    QgsProcessingContext,
    QgsProcessingFeedback,
    QgsRasterLayer,
    QgsRectangle,
    QgsCoordinateReferenceSystem,
)
import processing


class CanopyMapperAlgorithm(QgsProcessingAlgorithm):

    # --- Parametri: input raster ---
    ORTOFOTO = 'ORTOFOTO'
    RED_BAND = 'RED_BAND'
    GREEN_BAND = 'GREEN_BAND'
    BLUE_BAND = 'BLUE_BAND'
    HAS_NIR = 'HAS_NIR'
    NIR_BAND = 'NIR_BAND'

    # --- Parametri: soglia ---
    THRESHOLD_MODE = 'THRESHOLD_MODE'
    MANUAL_THRESHOLD = 'MANUAL_THRESHOLD'

    # --- Parametri: filtro poligoni ---
    MIN_AREA = 'MIN_AREA'
    MIN_COMPACTNESS = 'MIN_COMPACTNESS'

    # --- Parametri: altezza (opzionale) ---
    DSM = 'DSM'
    DTM = 'DTM'
    MIN_HEIGHT = 'MIN_HEIGHT'

    # --- Parametri: separazione chiome adiacenti (watershed) ---
    ENABLE_WATERSHED = 'ENABLE_WATERSHED'
    MIN_CROWN_DISTANCE = 'MIN_CROWN_DISTANCE'

    # --- Output ---
    OUTPUT = 'OUTPUT'

    def tr(self, string):
        return QCoreApplication.translate('CanopyMapper', string)

    def createInstance(self):
        return CanopyMapperAlgorithm()

    def name(self):
        return 'canopy_mapper'

    def displayName(self):
        return self.tr('Individua canopy da ortofoto (CanopyMapper)')

    def group(self):
        return self.tr('CanopyMapper')

    def groupId(self):
        return 'canopymapper'

    def shortHelpString(self):
        return self.tr(
            "Individua automaticamente le chiome di alberi e arbusti a "
            "partire da un'ortofoto (RGB o RGB+NIR), producendo poligoni "
            "vettoriali con area, perimetro, compattezza, classe "
            "dimensionale e, se fornisci DSM e DTM, altezza media/massima "
            "di ogni chioma (dal CHM = DSM - DTM).\n\n"
            "Se l'ortofoto ha una banda NIR usa l'NDVI (più affidabile); "
            "altrimenti calcola l'ExG (Excess Green) dalle sole bande "
            "RGB.\n\n"
            "Suggerimento: se non sai quale soglia usare, lascia la "
            "modalità 'Automatica (Otsu)': il metodo trova da solo il "
            "punto di separazione ottimale tra vegetazione e non.\n\n"
            "ATTENZIONE sulla soglia manuale: l'NDVI varia tra -1 e 1, "
            "l'ExG tipicamente tra -0.5 e 0.5. Un valore fuori da questo "
            "range (es. 2.0) non verrà mai superato da nessun pixel e "
            "produrrà un layer vuoto: non è un errore dello strumento.\n\n"
            "IMPORTANTE - alberi vicino a campi/prati: un indice di "
            "vegetazione da solo non distingue una chioma da un campo "
            "coltivato (sono verdi entrambi). Se un albero tocca un'area "
            "verde più ampia, i due vengono fusi in un unico poligono in "
            "fase di poligonizzazione, che poi il filtro di compattezza "
            "scarta insieme al campo. Se fornisci DSM e DTM e imposti "
            "un'altezza minima, questo problema si riduce molto, perché "
            "campi e prati vengono esclusi già in base alla quota, prima "
            "ancora della poligonizzazione.\n\n"
            "SEPARAZIONE CHIOME ADIACENTI (watershed): quando più chiome si "
            "toccano (es. un gruppo di alberi, o alberi accanto a una "
            "siepe), di norma verrebbero fuse in un unico poligono "
            "irregolare e scartate dal filtro di compattezza. Con "
            "'Separa chiome adiacenti' attivo, lo strumento individua le "
            "cime (massimi locali di altezza, se hai fornito DSM/DTM, "
            "altrimenti la forma stessa della macchia) e divide ogni "
            "gruppo in poligoni singoli, uno per cima, prima di applicare "
            "i filtri. Richiede la libreria scipy (inclusa di default in "
            "QGIS 3.x recenti); se non disponibile, questo passaggio viene "
            "saltato automaticamente con un avviso, il resto funziona "
            "comunque."
        )

    # ------------------------------------------------------------------
    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.ORTOFOTO, self.tr('Ortofoto (RGB o RGB+NIR)')))

        self.addParameter(QgsProcessingParameterBand(
            self.RED_BAND, self.tr('Banda Rosso'),
            parentLayerParameterName=self.ORTOFOTO, defaultValue=1))
        self.addParameter(QgsProcessingParameterBand(
            self.GREEN_BAND, self.tr('Banda Verde'),
            parentLayerParameterName=self.ORTOFOTO, defaultValue=2))
        self.addParameter(QgsProcessingParameterBand(
            self.BLUE_BAND, self.tr('Banda Blu'),
            parentLayerParameterName=self.ORTOFOTO, defaultValue=3))

        self.addParameter(QgsProcessingParameterBoolean(
            self.HAS_NIR,
            self.tr("L'ortofoto ha una banda infrarosso vicino (NIR) -> "
                    "usa NDVI invece di ExG"),
            defaultValue=False))
        nir_param = QgsProcessingParameterBand(
            self.NIR_BAND, self.tr('Banda NIR (se presente)'),
            parentLayerParameterName=self.ORTOFOTO, defaultValue=4,
            optional=True)
        self.addParameter(nir_param)

        self.addParameter(QgsProcessingParameterEnum(
            self.THRESHOLD_MODE,
            self.tr('Modalità soglia vegetazione'),
            options=[self.tr('Automatica (Otsu)'),
                     self.tr('Manuale')],
            defaultValue=0))
        self.addParameter(QgsProcessingParameterNumber(
            self.MANUAL_THRESHOLD,
            self.tr('Soglia manuale (usata solo se modalità = Manuale)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.2, optional=True))

        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_AREA,
            self.tr('Area minima poligono (m²) - scarta il rumore'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=1.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_COMPACTNESS,
            self.tr('Compattezza minima (0-1, 1=cerchio perfetto) - '
                     'scarta ombre/forme allungate innaturali'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.15, minValue=0.0, maxValue=1.0))

        dsm_param = QgsProcessingParameterRasterLayer(
            self.DSM, self.tr('DSM - Modello Superficie (opzionale, per altezza)'),
            optional=True)
        dtm_param = QgsProcessingParameterRasterLayer(
            self.DTM, self.tr('DTM - Modello Terreno (opzionale, per altezza)'),
            optional=True)
        self.addParameter(dsm_param)
        self.addParameter(dtm_param)
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_HEIGHT,
            self.tr('Altezza minima chioma (m) - richiede DSM e DTM. '
                     'Esclude a monte prati/campi bassi usando il CHM, '
                     'prima ancora della poligonizzazione (consigliato '
                     'per aree agricole/rurali)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=2.0, optional=True, minValue=0.0))

        self.addParameter(QgsProcessingParameterBoolean(
            self.ENABLE_WATERSHED,
            self.tr('Separa chiome adiacenti (watershed) - divide gruppi di '
                     'alberi fusi insieme in poligoni singoli invece di '
                     'scartarli in blocco'),
            defaultValue=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.MIN_CROWN_DISTANCE,
            self.tr('Distanza minima tra cime (m) - due massimi più vicini '
                     'di così vengono considerati la stessa chioma. Include '
                     'anche un leggero appianamento della superficie per '
                     'non scambiare micro-rumore del CHM per cime separate '
                     '(evita di frammentare una chioma singola in troppi '
                     'pezzi minuscoli, difficili da selezionare in QGIS)'),
            type=QgsProcessingParameterNumber.Double,
            defaultValue=4.0, minValue=0.5))

        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT, self.tr('Canopy risultante (poligoni)')))

    # ------------------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):

        ortofoto = self.parameterAsRasterLayer(parameters, self.ORTOFOTO, context)
        if ortofoto is None:
            raise QgsProcessingException(self.tr("Ortofoto non valida."))

        red_b = self.parameterAsInt(parameters, self.RED_BAND, context)
        green_b = self.parameterAsInt(parameters, self.GREEN_BAND, context)
        blue_b = self.parameterAsInt(parameters, self.BLUE_BAND, context)
        has_nir = self.parameterAsBoolean(parameters, self.HAS_NIR, context)
        nir_b = self.parameterAsInt(parameters, self.NIR_BAND, context) if has_nir else None

        threshold_mode = self.parameterAsEnum(parameters, self.THRESHOLD_MODE, context)
        manual_threshold = self.parameterAsDouble(parameters, self.MANUAL_THRESHOLD, context)

        min_area = self.parameterAsDouble(parameters, self.MIN_AREA, context)
        min_compactness = self.parameterAsDouble(parameters, self.MIN_COMPACTNESS, context)

        dsm = self.parameterAsRasterLayer(parameters, self.DSM, context)
        dtm = self.parameterAsRasterLayer(parameters, self.DTM, context)
        compute_height = dsm is not None and dtm is not None
        min_height = self.parameterAsDouble(parameters, self.MIN_HEIGHT, context) \
            if parameters.get(self.MIN_HEIGHT) is not None else None
        use_height_prefilter = compute_height and min_height is not None

        feedback.pushInfo(self.tr('--- STEP 1/7: calcolo indice di vegetazione ---'))
        index_path = self._compute_vegetation_index(
            ortofoto, red_b, green_b, blue_b, nir_b, has_nir, context, feedback)

        feedback.pushInfo(self.tr('--- STEP 2/7: sogliatura ---'))
        if threshold_mode == 0:
            threshold = self._otsu_threshold(index_path, feedback)
            feedback.pushInfo(self.tr(
                'Soglia automatica (Otsu) calcolata: {0:.4f}').format(threshold))
        else:
            threshold = manual_threshold
            feedback.pushInfo(self.tr('Soglia manuale usata: {0:.4f}').format(threshold))

        mask_path = self._threshold_to_mask(index_path, threshold, context, feedback)

        chm_path = None
        chm_aligned_path = None
        if compute_height:
            feedback.pushInfo(self.tr('--- STEP 3/7: calcolo CHM (DSM - DTM) ---'))
            chm_path = self._compute_chm(dsm, dtm, context, feedback)
            chm_aligned_path = self._align_raster_to_ortofoto(
                chm_path, ortofoto, context, feedback)

        if use_height_prefilter:
            feedback.pushInfo(self.tr(
                '--- STEP 3b/7: esclusione vegetazione bassa (< {0} m) prima della '
                'poligonizzazione ---').format(min_height))
            mask_path = self._apply_height_prefilter(
                mask_path, chm_aligned_path, min_height, context, feedback)
        else:
            feedback.pushInfo(self.tr(
                '--- STEP 3b/7: nessun filtro altezza a monte (DSM/DTM o altezza minima '
                'non impostati) ---'))

        feedback.pushInfo(self.tr('--- STEP 4/7: poligonizzazione ---'))
        is_labeled = False
        if self.parameterAsBoolean(parameters, self.ENABLE_WATERSHED, context):
            min_crown_dist = self.parameterAsDouble(parameters, self.MIN_CROWN_DISTANCE, context)
            pixel_size = ortofoto.rasterUnitsPerPixelX()
            watershed_path = self._watershed_split(
                mask_path, chm_aligned_path,
                min_crown_dist, pixel_size, context, feedback)
            if watershed_path is not None:
                mask_path = watershed_path
                is_labeled = True
        raw_polygons = self._polygonize(mask_path, is_labeled, context, feedback)

        feedback.pushInfo(self.tr('--- STEP 5/7: pulizia geometrie e calcolo area/perimetro/compattezza ---'))
        clean_polygons = self._clean_and_measure(raw_polygons, context, feedback)

        feedback.pushInfo(self.tr('--- STEP 6/7: filtro per area minima e compattezza ---'))
        filtered = self._filter_polygons(clean_polygons, min_area, min_compactness, context, feedback)

        if compute_height:
            feedback.pushInfo(self.tr('--- STEP 7/7: statistiche altezza per chioma (dal CHM) ---'))
            filtered = self._add_height_from_chm(filtered, chm_path, context, feedback)
        else:
            feedback.pushInfo(self.tr(
                '--- STEP 7/7: DSM/DTM non forniti, salto il calcolo altezza ---'))

        final_layer = self._finalize_attributes(filtered, compute_height, context, feedback)

        (sink, dest_id) = self.parameterAsSink(
            parameters, self.OUTPUT, context,
            final_layer.fields(), final_layer.wkbType(), final_layer.sourceCrs())
        if sink is None:
            raise QgsProcessingException(
                self.invalidSinkError(parameters, self.OUTPUT))

        for feat in final_layer.getFeatures():
            sink.addFeature(feat)

        feedback.pushInfo(self.tr('Completato. Chiome individuate: {0}').format(
            final_layer.featureCount()))

        return {self.OUTPUT: dest_id}

    # ------------------------------------------------------------------
    # STEP 1 - indice di vegetazione
    # ------------------------------------------------------------------
    def _compute_vegetation_index(self, ortofoto, red_b, green_b, blue_b,
                                   nir_b, has_nir, context, feedback):
        source = ortofoto.source()
        if has_nir and nir_b:
            # NDVI = (NIR - RED) / (NIR + RED)
            expression = (
                f'(A.astype(float) - B.astype(float)) / '
                f'((A.astype(float) + B.astype(float)) + 1e-6)'
            )
            feedback.pushInfo(self.tr('Uso NDVI (bande NIR/Rosso).'))
            result = processing.run('gdal:rastercalculator', {
                'INPUT_A': source, 'BAND_A': nir_b,
                'INPUT_B': source, 'BAND_B': red_b,
                'FORMULA': '(A.astype(float)-B.astype(float))/((A.astype(float)+B.astype(float))+0.000001)',
                'NO_DATA': None, 'RTYPE': 5,  # Float32
                'OUTPUT': 'TEMPORARY_OUTPUT',
            }, context=context, feedback=feedback)
        else:
            # ExG = 2*G - R - B, normalizzato in 0-1
            feedback.pushInfo(self.tr(
                'Nessuna banda NIR indicata: uso ExG (Excess Green) da RGB.'))
            result = processing.run('gdal:rastercalculator', {
                'INPUT_A': source, 'BAND_A': green_b,
                'INPUT_B': source, 'BAND_B': red_b,
                'INPUT_C': source, 'BAND_C': blue_b,
                'FORMULA': '(2*A.astype(float)-B.astype(float)-C.astype(float))/'
                           '(A.astype(float)+B.astype(float)+C.astype(float)+0.000001)',
                'NO_DATA': None, 'RTYPE': 5,
                'OUTPUT': 'TEMPORARY_OUTPUT',
            }, context=context, feedback=feedback)
        return result['OUTPUT']

    # ------------------------------------------------------------------
    # STEP 2a - soglia automatica di Otsu
    # ------------------------------------------------------------------
    def _otsu_threshold(self, raster_path, feedback):
        layer = QgsRasterLayer(raster_path, 'index')
        provider = layer.dataProvider()
        extent = layer.extent()
        width, height = layer.width(), layer.height()

        # Limitiamo la lettura per non esaurire la memoria su ortofoto enormi
        max_dim = 2000
        scale = min(1.0, max_dim / max(width, height))
        sample_w = max(1, int(width * scale))
        sample_h = max(1, int(height * scale))

        block = provider.block(1, extent, sample_w, sample_h)
        data = np.frombuffer(block.data(), dtype=np.float32).reshape(sample_h, sample_w).copy()
        nodata = provider.sourceNoDataValue(1)
        valid = data[np.isfinite(data)]
        if nodata is not None:
            valid = valid[valid != nodata]
        valid = valid[(valid > -1.0) & (valid < 1.0)]

        if valid.size == 0:
            feedback.pushWarning(self.tr(
                'Impossibile calcolare Otsu automaticamente, uso soglia 0.2 di default.'))
            return 0.2

        hist, bin_edges = np.histogram(valid, bins=256)
        hist = hist.astype(float)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        weight1 = np.cumsum(hist)
        weight2 = np.cumsum(hist[::-1])[::-1]
        # evita divisioni per zero
        weight1_safe = np.where(weight1 == 0, 1, weight1)
        weight2_safe = np.where(weight2 == 0, 1, weight2)

        mean1 = np.cumsum(hist * bin_centers) / weight1_safe
        mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2_safe[::-1])[::-1]

        variance_between = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
        if variance_between.size == 0:
            return 0.2
        idx = np.argmax(variance_between)
        return float(bin_centers[idx])

    # ------------------------------------------------------------------
    # STEP 2b - sogliatura -> maschera binaria
    # ------------------------------------------------------------------
    def _threshold_to_mask(self, index_path, threshold, context, feedback):
        result = processing.run('gdal:rastercalculator', {
            'INPUT_A': index_path, 'BAND_A': 1,
            'FORMULA': f'(A > {threshold}) * 1',
            'NO_DATA': 0, 'RTYPE': 0,  # Byte
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)
        return result['OUTPUT']

    # ------------------------------------------------------------------
    # STEP 3 - poligonizzazione
    # is_labeled=False: maschera binaria 0/1, teniamo dove valore==1
    # is_labeled=True: raster con un ID intero unico per chioma (da
    #   watershed), teniamo tutto tranne lo sfondo (valore==0)
    # ------------------------------------------------------------------
    def _polygonize(self, mask_path, is_labeled, context, feedback):
        result = processing.run('gdal:polygonize', {
            'INPUT': mask_path, 'BAND': 1, 'FIELD': 'valore',
            'EIGHT_CONNECTEDNESS': True, 'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)
        expr = '"valore" != 0' if is_labeled else '"valore" = 1'
        veg_only = processing.run('native:extractbyexpression', {
            'INPUT': result['OUTPUT'], 'EXPRESSION': expr,
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)
        return veg_only['OUTPUT']

    # ------------------------------------------------------------------
    # STEP 4 - pulizia geometrie + area/perimetro/compattezza
    # ------------------------------------------------------------------
    def _clean_and_measure(self, layer, context, feedback):
        fixed = processing.run('native:fixgeometries', {
            'INPUT': layer, 'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        with_geom = processing.run('qgis:exportaddgeometrycolumns', {
            'INPUT': fixed, 'CALC_METHOD': 0,  # coordinate del layer (metriche se CRS proiettato)
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        # rinomina/aggiunge compattezza = 4*pi*area / perimetro^2
        with_compactness = processing.run('native:fieldcalculator', {
            'INPUT': with_geom,
            'FIELD_NAME': 'compattez',
            'FIELD_TYPE': 0,  # Float
            'FIELD_LENGTH': 10, 'FIELD_PRECISION': 4,
            'FORMULA': '4 * pi() * "area" / ("perimeter" ^ 2)',
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        return with_compactness

    # ------------------------------------------------------------------
    # STEP 5 - filtro area minima e compattezza minima
    # ------------------------------------------------------------------
    def _filter_polygons(self, layer, min_area, min_compactness, context, feedback):
        expr = f'"area" >= {min_area} AND "compattez" >= {min_compactness}'
        result = processing.run('native:extractbyexpression', {
            'INPUT': layer, 'EXPRESSION': expr, 'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)
        return result['OUTPUT']

    # ------------------------------------------------------------------
    # CHM (DSM - DTM). Richiede che DSM e DTM condividano la stessa griglia
    # (stesso pixel size/estensione): è il caso tipico quando entrambi
    # derivano dallo stesso rilievo LiDAR/fotogrammetrico.
    # ------------------------------------------------------------------
    def _compute_chm(self, dsm, dtm, context, feedback):
        chm = processing.run('gdal:rastercalculator', {
            'INPUT_A': dsm.source(), 'BAND_A': 1,
            'INPUT_B': dtm.source(), 'BAND_B': 1,
            'FORMULA': '(A.astype(float) - B.astype(float))',
            'NO_DATA': None, 'RTYPE': 5,
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']
        return chm

    # ------------------------------------------------------------------
    # Riallinea (warp) un raster (tipicamente il CHM) sulla griglia
    # dell'ortofoto: stessa estensione, stessa risoluzione, stesso CRS.
    # Necessario perché DSM/DTM in generale hanno risoluzione/estensione
    # diverse dall'ortofoto, e qui servono confronti pixel-per-pixel.
    # ------------------------------------------------------------------
    def _align_raster_to_ortofoto(self, raster_path, ortofoto, context, feedback):
        extent_str = '{0},{1},{2},{3} [{4}]'.format(
            ortofoto.extent().xMinimum(), ortofoto.extent().xMaximum(),
            ortofoto.extent().yMinimum(), ortofoto.extent().yMaximum(),
            ortofoto.crs().authid())
        aligned = processing.run('gdal:warpreproject', {
            'INPUT': raster_path,
            'SOURCE_CRS': None,
            'TARGET_CRS': ortofoto.crs(),
            'RESAMPLING': 1,  # bilinear
            'TARGET_RESOLUTION': ortofoto.rasterUnitsPerPixelX(),
            'NODATA': -9999,
            'TARGET_EXTENT': extent_str,
            'TARGET_EXTENT_CRS': ortofoto.crs(),
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']
        return aligned

    # ------------------------------------------------------------------
    # Pre-filtro altezza: esclude dalla maschera vegetazione i pixel sotto
    # l'altezza minima (es. prati, campi coltivati) PRIMA della
    # poligonizzazione, così non si fondono più con gli alberi adiacenti.
    # Richiede un CHM già riallineato sulla griglia dell'ortofoto.
    # ------------------------------------------------------------------
    def _apply_height_prefilter(self, mask_path, chm_aligned_path,
                                 min_height, context, feedback):
        combined = processing.run('gdal:rastercalculator', {
            'INPUT_A': mask_path, 'BAND_A': 1,
            'INPUT_B': chm_aligned_path, 'BAND_B': 1,
            'FORMULA': f'((A == 1) & (B > {min_height})) * 1',
            'NO_DATA': 0, 'RTYPE': 0,
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']
        return combined

    # ------------------------------------------------------------------
    # Separazione chiome adiacenti (watershed "fatto in casa", solo numpy
    # + scipy, niente scikit-image). Individua le cime (massimi locali
    # dell'altezza se disponibile, altrimenti della distanza dal bordo
    # della macchia) ed espande ogni cima verso l'esterno assegnando ad
    # ogni pixel l'etichetta della cima raggiunta con priorità (quota più
    # alta prima) - così un gruppo di chiome fuse diventa tanti poligoni
    # separati invece di uno solo da accettare/scartare in blocco.
    # Ritorna il path del raster etichettato, oppure None se scipy non è
    # disponibile o non vengono trovate cime (in quel caso si prosegue
    # con la maschera binaria originale, comportamento pre-watershed).
    # ------------------------------------------------------------------
    def _watershed_split(self, mask_path, chm_path, min_distance_m,
                          pixel_size, context, feedback):
        try:
            from scipy import ndimage
        except ImportError:
            feedback.pushWarning(self.tr(
                "scipy non disponibile: salto la separazione delle chiome "
                "adiacenti (watershed). Il resto dello strumento funziona "
                "comunque, ma i gruppi di chiome fuse verranno gestiti solo "
                "dal filtro di compattezza, come nella versione precedente."))
            return None

        try:
            from osgeo import gdal
        except ImportError:
            feedback.pushWarning(self.tr(
                "Impossibile importare osgeo.gdal: salto la separazione "
                "watershed."))
            return None

        mask_ds = gdal.Open(mask_path)
        mask = mask_ds.GetRasterBand(1).ReadAsArray().astype(bool)
        geotransform = mask_ds.GetGeoTransform()
        projection = mask_ds.GetProjection()
        rows, cols = mask.shape

        if not mask.any():
            mask_ds = None
            return None

        if chm_path is not None:
            chm_ds = gdal.Open(chm_path)
            surface = chm_ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
            chm_ds = None
            surface = np.where(mask, surface, -1e6)
        else:
            surface = ndimage.distance_transform_edt(mask).astype(np.float32)

        min_distance_px = max(1, int(round(min_distance_m / pixel_size)))
        smoothing_sigma = min_distance_px * 0.5
        if smoothing_sigma > 0:
            filled = np.where(mask, surface, 0)
            smoothed = ndimage.gaussian_filter(filled, sigma=smoothing_sigma)
            surface = np.where(mask, smoothed, -1e6).astype(np.float32)
        footprint = np.ones((2 * min_distance_px + 1, 2 * min_distance_px + 1))
        local_max = (ndimage.maximum_filter(surface, footprint=footprint) == surface) & mask
        markers, n_markers = ndimage.label(local_max, structure=np.ones((3, 3)))

        feedback.pushInfo(self.tr(
            'Watershed: {0} cime individuate come marker').format(n_markers))
        if n_markers == 0:
            mask_ds = None
            return None

        labels = self._priority_flood_watershed(surface, markers, mask)

        out_path = mask_path[:-4] + '_watershed.tif' if mask_path.lower().endswith('.tif') \
            else mask_path + '_watershed.tif'
        driver = gdal.GetDriverByName('GTiff')
        out_ds = driver.Create(out_path, cols, rows, 1, gdal.GDT_Int32)
        out_ds.SetGeoTransform(geotransform)
        out_ds.SetProjection(projection)
        out_band = out_ds.GetRasterBand(1)
        out_band.WriteArray(labels.astype(np.int32))
        out_band.SetNoDataValue(0)
        out_band.FlushCache()
        out_ds = None
        mask_ds = None
        return out_path

    # ------------------------------------------------------------------
    # Nucleo dell'algoritmo: espansione a priorità (flood) dai marker.
    # Ad ogni passo si assegna al vicino non ancora etichettato con la
    # quota più alta disponibile in coda; equivalente a un watershed
    # marker-controlled sulla superficie 'surface', limitato a 'mask'.
    # ------------------------------------------------------------------
    def _priority_flood_watershed(self, surface, markers, mask):
        rows, cols = surface.shape
        labels = markers.copy()
        heap = []
        counter = 0

        def push_neighbors(r, c, label):
            nonlocal counter
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols and mask[nr, nc] \
                            and labels[nr, nc] == 0:
                        counter += 1
                        heapq.heappush(heap, (-surface[nr, nc], counter, nr, nc, label))

        seed_rows, seed_cols = np.nonzero(markers)
        for r, c in zip(seed_rows, seed_cols):
            push_neighbors(r, c, markers[r, c])

        while heap:
            _, _, r, c, label = heapq.heappop(heap)
            if labels[r, c] != 0:
                continue
            labels[r, c] = label
            push_neighbors(r, c, label)

        return labels

    # ------------------------------------------------------------------
    # Statistiche zonali di altezza per chioma (usa il CHM alla sua
    # risoluzione/griglia originale: QGIS gestisce la corrispondenza
    # vettore/raster anche con CRS dichiarati diversi).
    # ------------------------------------------------------------------
    def _add_height_from_chm(self, polygons, chm_path, context, feedback):
        result = processing.run('native:zonalstatisticsfb', {
            'INPUT': polygons, 'INPUT_RASTER': chm_path, 'RASTER_BAND': 1,
            'COLUMN_PREFIX': 'h_',
            'STATISTICS': [2, 6],  # 2=media, 6=massimo
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)
        return result['OUTPUT']

    # ------------------------------------------------------------------
    # Arrotonda un campo numerico esistente a N decimali (aggiorna il
    # campo stesso, non ne crea uno nuovo).
    # ------------------------------------------------------------------
    def _round_field(self, layer, field_name, decimals, context, feedback):
        return processing.run('native:fieldcalculator', {
            'INPUT': layer,
            'FIELD_NAME': field_name,
            'FIELD_TYPE': 0,  # Float
            'FIELD_LENGTH': 12, 'FIELD_PRECISION': decimals,
            'FORMULA': f'round("{field_name}", {decimals})',
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

    # ------------------------------------------------------------------
    # Finalizzazione: arrotondamento campi, ID progressivo + classe
    # dimensionale automatica
    # ------------------------------------------------------------------
    def _finalize_attributes(self, layer, has_height, context, feedback):
        feedback.pushInfo(self.tr('Arrotondamento campi numerici a 2 decimali...'))
        rounded = layer
        for fname in ('area', 'perimeter', 'compattez'):
            rounded = self._round_field(rounded, fname, 2, context, feedback)
        if has_height:
            for fname in ('h_mean', 'h_max'):
                rounded = self._round_field(rounded, fname, 2, context, feedback)
        layer = rounded

        with_id = processing.run('native:addautoincrementalfield', {
            'INPUT': layer, 'FIELD_NAME': 'id_pianta', 'START': 1,
            'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        # Classificazione dimensionale automatica, ora incrociata con
        # l'altezza quando disponibile (non solo l'area). Senza questo
        # incrocio, un'area molto estesa ma bassa (es. un vigneto/frutteto
        # rimasto residualmente nella maschera, con chiome non fuse tra
        # loro ma sparse su una superficie ampia) verrebbe etichettata
        # "albero_grande" solo perché occupa molti m², anche se in realtà
        # è alta pochi metri - inverosimile per un albero grande vero.
        # Soglie di partenza generiche: tarale sul tuo contesto (specie,
        # età, urbano/extraurbano) editando i valori qui sotto.
        if has_height:
            formula = (
                "CASE "
                "WHEN \"area\" < 4 THEN 'arbusto' "
                "WHEN \"area\" < 15 THEN 'albero_piccolo' "
                "WHEN \"area\" < 40 AND \"h_mean\" >= 3 THEN 'albero_medio' "
                "WHEN \"area\" >= 40 AND \"h_mean\" >= 6 THEN 'albero_grande' "
                "WHEN \"h_mean\" < 4 THEN 'sospetto_area_ampia_bassa' "
                "ELSE 'albero_medio' END"
            )
        else:
            formula = (
                "CASE "
                "WHEN \"area\" < 4 THEN 'arbusto' "
                "WHEN \"area\" < 15 THEN 'albero_piccolo' "
                "WHEN \"area\" < 40 THEN 'albero_medio' "
                "ELSE 'albero_grande' END"
            )
        classified = processing.run('native:fieldcalculator', {
            'INPUT': with_id, 'FIELD_NAME': 'classe_dim',
            'FIELD_TYPE': 2, 'FIELD_LENGTH': 30, 'FIELD_PRECISION': 0,
            'FORMULA': formula, 'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        # Campo vuoto pronto per la classificazione manuale della specie
        with_species = processing.run('native:fieldcalculator', {
            'INPUT': classified, 'FIELD_NAME': 'specie',
            'FIELD_TYPE': 2, 'FIELD_LENGTH': 50, 'FIELD_PRECISION': 0,
            'FORMULA': "''", 'OUTPUT': 'TEMPORARY_OUTPUT',
        }, context=context, feedback=feedback)['OUTPUT']

        return with_species
