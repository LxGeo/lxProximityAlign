#include"gdal_algs/proximity_transformer.h"
#include <gdal_priv.h>
#include "gdal_utils.h"
#include "gdal_alg.h"
#include "io_raster.h"


namespace LxGeo
{
	namespace lxProximityAlign
	{
		using namespace IO_DATA;

		void transform_to_proximity(const std::string& in_raster_path, const std::string& out_raster_path) {

			RasterIO in_raster = RasterIO(in_raster_path, GA_ReadOnly, true);
			GDALRasterBandH in_band = in_raster.raster_dataset->GetRasterBand(1);

			GDALDataset* out_dataset = in_raster.create_copy_dataset(out_raster_path, GDT_Float32, 1);
			GDALRasterBandH out_band = out_dataset->GetRasterBand(1);

			char** argv = NULL;

			argv = CSLAddString(argv, "-DISTUNITS");
			argv = CSLAddString(argv, "PIXEL");

			argv = CSLAddString(argv, "-maxdist");
			argv = CSLAddString(argv, "500");

			argv = CSLAddString(argv, "-nodata");
			argv = CSLAddString(argv, "0");

			GDALComputeProximity(in_band, out_band, argv, NULL, NULL);

			CSLDestroy(argv);
			GDALClose(out_dataset);

		}

	}
}