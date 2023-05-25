#include "gdal_algs/rasterizer.h"


namespace LxGeo
{
    namespace lxProximityAlign
    {

        void rasterize_shapefile(const std::string& out_raster_path, const std::string& input_shapefile_path, OGREnvelope* raster_extents, double raster_px_size, double raster_py_size)
        {
            unsigned int openFlags = GDAL_OF_VECTOR | GDAL_OF_READONLY;
            GDALDataset* pSrcDataset = static_cast<GDALDataset*>(GDALOpenEx(input_shapefile_path.c_str(), openFlags, NULL, NULL, NULL));
            OGRLayer* pSrcLayer = pSrcDataset->GetLayer(0);

            // Rasterize options.
            char** argv = NULL;

            argv = CSLAddString(argv, "-burn");
            argv = CSLAddString(argv, "1");

            argv = CSLAddString(argv, "-l");
            argv = CSLAddString(argv, pSrcLayer->GetName());

            argv = CSLAddString(argv, "-a_nodata");
            argv = CSLAddString(argv, "-32600");

            argv = CSLAddString(argv, "-tr");
            argv = CSLAddString(argv, std::to_string(raster_px_size).c_str());
            argv = CSLAddString(argv, std::to_string(raster_py_size).c_str());

            argv = CSLAddString(argv, "-ot");
            argv = CSLAddString(argv, "Byte");

            argv = CSLAddString(argv, "-at");

            if (raster_extents) {
                argv = CSLAddString(argv, "-te");
                argv = CSLAddString(argv, std::to_string(raster_extents->MinX-1).c_str());
                argv = CSLAddString(argv, std::to_string(raster_extents->MinY-1).c_str());
                argv = CSLAddString(argv, std::to_string(raster_extents->MaxX+1).c_str());
                argv = CSLAddString(argv, std::to_string(raster_extents->MaxY+1).c_str());
            }
            else {
                OGREnvelope out_extents; pSrcLayer->GetExtent(&out_extents);
                argv = CSLAddString(argv, "-te");
                argv = CSLAddString(argv, std::to_string(out_extents.MinX).c_str());
                argv = CSLAddString(argv, std::to_string(out_extents.MinY).c_str());
                argv = CSLAddString(argv, std::to_string(out_extents.MaxX).c_str());
                argv = CSLAddString(argv, std::to_string(out_extents.MaxY).c_str());
            }

            GDALRasterizeOptions* pOptions = GDALRasterizeOptionsNew(argv, NULL);

            // Perform rasterization.
            int          usageError;
            GDALDataset* pDstDataset =
                static_cast<GDALDataset*>(GDALRasterize(out_raster_path.c_str(), NULL, pSrcDataset, pOptions, &usageError));

            //------------------------------
            // Cleanup.
            GDALRasterizeOptionsFree(pOptions);
            CSLDestroy(argv);
            GDALClose(pSrcDataset);
            GDALClose(pDstDataset);
        }

    }
}