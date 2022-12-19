context ("data_frame-osm")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

test_that ("multipolygon", {
    osm_multi <- test_path ("fixtures", "osm-multi.osm")
    x_sf <- sf::st_read (
        osm_multi,
        layer = "multipolygons",
        stringsAsFactors = FALSE,
        quiet = TRUE
    )
    q0 <- opq (bbox = c (1, 1, 5, 5))
    x <- osmdata_data_frame (q0, osm_multi)

    # GDAL spits out a whole lot of generic field names, so first the
    # two have to be reduced to common fields.
    x_sf <- sf::st_drop_geometry (x_sf)
    x <- x [, which (names (x) %in% names (x_sf))]
    x_sf <- x_sf [, which (names (x_sf) %in% names (x))]

    expect_identical (names (x), names (x_sf))
    expect_true (x_sf$osm_id %in% x$osm_id)
    expect_true (x_sf$name %in% x$name)
    expect_true (x_sf$type %in% x$type)
})


test_that ("multilinestring", {
    osm_multi <- test_path ("fixtures", "osm-multi.osm")
    x_sf <- sf::st_read (
        osm_multi,
        layer = "multilinestrings",
        stringsAsFactors = FALSE,
        quiet = TRUE
    )
    q0 <- opq (bbox = c (1, 1, 5, 5))
    x <- osmdata_data_frame (q0, osm_multi)
    x_sf <- sf::st_drop_geometry (x_sf)
    x <- x [, which (names (x) %in% names (x_sf))]
    x_sf <- x_sf [, which (names (x_sf) %in% names (x))]

    expect_identical (names (x), names (x_sf))
    expect_true (x_sf$osm_id %in% x$osm_id)
    expect_true (x_sf$name %in% x$name)
    expect_true (x_sf$type %in% x$type)
})

test_that ("ways", {
    osm_ways <- test_path ("fixtures", "osm-ways.osm")
    x_sf <- sf::st_read (
        osm_ways,
        layer = "lines",
        stringsAsFactors = FALSE,
        quiet = TRUE
    )
    q0 <- opq (bbox = c (1, 1, 5, 5))
    x <- osmdata_data_frame (q0, osm_ways)
    x_sf <- sf::st_drop_geometry (x_sf)
    x <- x [, which (names (x) %in% names (x_sf))]
    x_sf <- x_sf [, which (names (x_sf) %in% names (x))]

    expect_setequal (names (x), names (x_sf))
    expect_true (all (x_sf$osm_id %in% x$osm_id))
    expect_true (all (x_sf$name %in% x$name))
    expect_true (all (x_sf$highway %in% x$highway))
})

test_that ("empty result", {
    q0 <- getbb ("Països Catalans", featuretype = "relation") %>%
        opq (nodes_only = TRUE, datetime = "1714-09-11T00:00:00Z") %>%
        add_osm_feature ("does not exist", "&%$")

    osm_empty <- test_path ("fixtures", "osm-empty.osm")
    doc <- xml2::read_xml (osm_empty)

    x <- osmdata_data_frame (q0, doc)

    cols <- c ("osm_type", "osm_id")
    expect_named (x, cols)
    expect_s3_class (x, "data.frame")
    expect_identical (nrow (x), 0L)

    obj_overpass_call <- osmdata (bbox = q0$bbox, overpass_call = opq_string_intern (q0))
    obj_opq <- osmdata (bbox = q0$bbox, overpass_call = q0)
    obj <- osmdata (bbox = q0$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = get_metadata (obj, doc)$meta)

    expect_equal (metaL$meta_overpass_call$query_type, "date")
    expect_equal (metaL$meta_opq$query_type, "date")
    expect_null (metaL$meta_no_call$query_type)

    # adiff
    q0 <- getbb ("Països Catalans", featuretype = "relation") %>%
        opq (nodes_only = TRUE, datetime = "1714-09-11T00:00:00Z", adiff = TRUE) %>%
        add_osm_feature ("does not exist", "&%$")

    # osm_empty <- test_path ("fixtures", "osm-empty.osm") # same result
    # doc <- xml2::read_xml (osm_empty)

    x <- osmdata_data_frame (q0, doc)

    cols <- c ("osm_type", "osm_id")
    expect_named (x, cols)
    expect_s3_class (x, "data.frame")
    expect_identical (nrow (x), 0L)

    obj_overpass_call <- osmdata (bbox = q0$bbox, overpass_call = opq_string_intern (q0))
    obj_opq <- osmdata (bbox = q0$bbox, overpass_call = q0)
    obj <- osmdata (bbox = q0$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = get_metadata (obj, doc)$meta)

    expect_equal (metaL$meta_overpass_call$query_type, "adiff")
    expect_equal (metaL$meta_opq$query_type, "adiff")
    expect_null (metaL$meta_no_call$query_type)

    expect_identical (metaL$meta_overpass_call, metaL$meta_opq)
    expect_identical (metaL$meta_overpass_call$datetime_from, attr (q0, "datetime"))
})

test_that ("attributes", {
    osm_multi <- test_path ("fixtures", "osm-multi.osm")

    q0 <- opq (bbox = c (1, 1, 5, 5))
    x <- osmdata_data_frame (q0, osm_multi)
    x_no_call <- osmdata_data_frame (doc = osm_multi)
    x_sf <- osmdata_sf (q0, osm_multi)

    expect_s3_class (x, "data.frame")
    expect_identical (names (attributes (x)),
                      c ("names", "class", "row.names", "bbox", "overpass_call", "meta"))
    expect_identical (attr (x, "bbox"), q0$bbox)
    expect_identical (attr (x, "overpass_call"), x_sf$overpass_call)
    expect_identical (attr (x, "meta"), x_sf$meta)
    # no call
    expect_s3_class (x_no_call, "data.frame")
    expect_identical (names (attributes (x_no_call)),
                      c ("names", "class", "row.names", "meta"))
    expect_identical (attr (x_no_call, "meta"), x_sf$meta)
})

test_that ("date", {
    q <- getbb ("Conflent", featuretype = "relation") %>%
        opq (nodes_only = TRUE, datetime = "2020-11-07T00:00:00Z") %>%
        add_osm_feature ("natural", "peak") %>%
        add_osm_feature ("prominence")  %>%
        add_osm_feature ("name:ca")

    osm_meta_date <- test_path ("fixtures", "osm-date.osm")
    doc <- xml2::read_xml (osm_meta_date)

    x <- osmdata_data_frame (q, doc)
    x_no_call <- osmdata_data_frame (doc = doc, quiet = FALSE)

    cols <- c ("osm_type", "osm_id", "ele", "name",
               "name:ca", "natural", "prominence")
    expect_named (x, cols)
    expect_named (x_no_call, cols)
    expect_s3_class (x, "data.frame")
    expect_s3_class (x_no_call, "data.frame")

    obj_overpass_call <- osmdata (bbox = q$bbox, overpass_call = opq_string_intern (q))
    obj_opq <- osmdata (bbox = q$bbox, overpass_call = q)
    obj <- osmdata (bbox = q$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = get_metadata (obj, doc)$meta)

    expect_identical (metaL$meta_overpass_call, metaL$meta_opq)
    expect_identical (metaL$meta_overpass_call$datetime_to, attr (q, "datetime"))
    expect_null (metaL$meta_overpass_call$datetime_from)
    expect_null (metaL$meta_no_call$query_type)
})

test_that ("out meta & diff", {
    q <- getbb ("Conflent", featuretype = "relation") %>%
        opq (nodes_only = TRUE, out = "meta", datetime = "2020-11-07T00:00:00Z",
             datetime2 = "2022-12-04T00:00:00Z") %>%
        add_osm_feature ("natural", "peak") %>%
        add_osm_feature ("prominence")  %>%
        add_osm_feature ("name:ca")

    osm_meta_diff <- test_path ("fixtures", "osm-meta_diff.osm")
    doc <- xml2::read_xml (osm_meta_diff)

    x <- osmdata_data_frame (q, doc, quiet = FALSE)
    x_no_call <- osmdata_data_frame (doc = doc)

    cols <- c ("osm_type", "osm_id", "osm_version", "osm_timestamp",
               "osm_changeset", "osm_uid", "osm_user", "ele", "name",
               "name:ca", "natural", "prominence")
    expect_named (x, cols)
    expect_named (x_no_call, cols)
    expect_s3_class (x, "data.frame")
    expect_s3_class (x_no_call, "data.frame")

    obj_overpass_call <- osmdata (bbox = q$bbox, overpass_call = opq_string_intern (q))
    obj_opq <- osmdata (bbox = q$bbox, overpass_call = q)
    obj <- osmdata (bbox = q$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = get_metadata (obj, doc)$meta)

    expect_identical (metaL$meta_overpass_call, metaL$meta_opq)
    expect_identical (metaL$meta_overpass_call$datetime_from, attr (q, "datetime"))
    expect_identical (metaL$meta_overpass_call$datetime_to, attr (q, "datetime2"))
    expect_null (metaL$meta_no_call$query_type)
})

test_that ("out meta & adiff", {
    q <- getbb ("Conflent", featuretype = "relation") %>%
        opq (nodes_only = TRUE, out = "meta",
             datetime = "2020-11-07T00:00:00Z", adiff = TRUE) %>%
        add_osm_feature ("natural", "peak") %>%
        add_osm_feature ("prominence")  %>%
        add_osm_feature ("name:ca")

    osm_meta_adiff <- test_path ("fixtures", "osm-meta_adiff.osm")
    doc <- xml2::read_xml (osm_meta_adiff)

    expect_silent( x <- osmdata_data_frame (opq_string_intern (q), doc) )
    expect_warning(
        x_no_call <- osmdata_data_frame (doc = doc),
        "OSM data is ambiguous and can correspond either to a diff or an adiff query."
    ) # query_type assigned to diff

    cols <- c ("osm_type", "osm_id",
               "osm_version", "osm_timestamp", "osm_changeset", "osm_uid", "osm_user",
               "adiff_action", "adiff_date", "adiff_visible",
               "ele", "name", "name:ca", "natural", "prominence",
               "source:prominence", "wikidata", "wikipedia")
    expect_named (x, cols)
    # expect_named (x_no_call, cols) # query_type assigned to diff
    expect_s3_class (x, "data.frame")
    expect_s3_class (x_no_call, "data.frame")

    obj_overpass_call <- osmdata (bbox = q$bbox, overpass_call = opq_string_intern (q))
    obj_opq <- osmdata (bbox = q$bbox, overpass_call = q)
    obj <- osmdata (bbox = q$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = expect_warning (
                                     get_metadata (obj, doc)$meta,
                                     "OSM data is ambiguous and can correspond either to a diff or an adiff query"
                                 )
                  )

    expect_equal (metaL$meta_overpass_call$query_type, "adiff")
    expect_equal (metaL$meta_opq$query_type, "adiff")
    expect_equal (metaL$meta_no_call$query_type, "diff")
    expect_identical (metaL$meta_overpass_call, metaL$meta_opq)
    expect_identical (metaL$meta_overpass_call$datetime_from, attr (q, "datetime"))
})

test_that ("adiff2", {
    q <- getbb ("Perpinyà", featuretype = "relation") %>%
        opq (nodes_only = TRUE,
             datetime = "2012-11-07T00:00:00Z",
             datetime2 = "2016-11-07T00:00:00Z",
             adiff = TRUE) %>%
        add_osm_feature ("amenity", "restaurant")

    osm_adiff2 <- test_path ("fixtures", "osm-adiff2.osm")
    doc <- xml2::read_xml (osm_adiff2)

    x <- osmdata_data_frame (q, doc, quiet = FALSE)
    x_no_call <- osmdata_data_frame (doc = doc)

    cols <- c ("osm_type", "osm_id",
               "adiff_action", "adiff_date", "adiff_visible",
               "addr:housenumber", "addr:street", "amenity", "created_by",
               "cuisine", "name", "phone")
    expect_named (x, cols)
    expect_named (x_no_call, cols)
    expect_s3_class (x, "data.frame")
    expect_s3_class (x_no_call, "data.frame")

    obj_overpass_call <- osmdata (bbox = q$bbox, overpass_call = opq_string_intern (q))
    obj_opq <- osmdata (bbox = q$bbox, overpass_call = q)
    obj <- osmdata (bbox = q$bbox)

    metaL<- list (meta_overpass_call = get_metadata (obj_overpass_call, doc)$meta,
                  meta_opq = get_metadata (obj_opq, doc)$meta,
                  meta_no_call = get_metadata (obj, doc)$meta)

    k <- sapply (metaL, function (x) expect_equal (x$query_type, "adiff"))
    expect_identical (metaL$meta_overpass_call, metaL$meta_opq)
    expect_identical (metaL$meta_overpass_call$datetime_from, attr (q, "datetime"))
})