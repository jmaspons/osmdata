/***************************************************************************
 *  Project:    osmdata
 *  File:       osmdata-df.cpp
 *  Language:   C++
 *
 *  osmdata is free software: you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free
 *  Software Foundation, either version 3 of the License, or (at your option)
 *  any later version.
 *
 *  osmdata is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  osm-router.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Extract OSM data from an object of class XmlData and return
 *                  it in Rcpp::List format.
 *
 *  Limitations:
 *
 *  Dependencies:       none (rapidXML header included in osmdata)
 *
 *  Compiler Options:   -std=c++11
 ***************************************************************************/

#include "osmdata.h"

#include <Rcpp.h>

// Note: This code uses explicit index counters within most loops which use Rcpp
// objects, because these otherwise require a
// static_cast <size_t> (std::distance (...)). This operation copies each
// instance and can slow the loops down by several orders of magnitude!

/************************************************************************
 ************************************************************************
 **                                                                    **
 **          1. PRIMARY FUNCTIONS TO TRACE WAYS AND RELATIONS          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/


//' get_osm_relations
//'
//' Return a dual Rcpp::List containing all OSM relations tags.
//'
//' @param rels Pointer to the vector of Relation objects
//' @param unique_vals Pointer to a UniqueVals object containing std::sets of all
//'       unique IDs and keys for each kind of OSM object (nodes, ways, rels).
//'
//' @return A dual Rcpp::List, the first of which contains the multipolygon
//'         relations; the second the multilinestring relations.
//'
//' @noRd
Rcpp::List osm_df::get_osm_relations (const Relations &rels,
        const UniqueVals &unique_vals,
        const Rcpp::NumericVector &bbox)
{
    /* Trace all multipolygon relations. These are the only OSM types where
     * sizes are not known before, so lat-lons and node names are stored in
     * dynamic vectors. These are 3D monsters: #1 for relation, #2 for polygon
     * in relation, and #3 for data. There are also associated 2D vector<vector>
     * objects for IDs and multilinestring roles. */
    std::set <std::string> keyset; // must be ordered!
    std::vector <std::string> colnames = {"lat", "lon"}, rownames;
    Rcpp::List dimnames (0);
    Rcpp::NumericMatrix nmat (Rcpp::Dimension (0, 0));

    double_arr2 lat_vec, lon_vec;
    double_arr3 lat_arr_rel, lon_arr_rel;
    string_arr2 rowname_vec, id_vec_rel;
    string_arr3 rowname_arr_mp;
    std::vector <osmid_t> ids_ls;
    std::vector <std::string> ids_rel, rel_id;
    osmt_arr2 id_vec_ls;
    std::vector <std::string> roles;

    unsigned int nmp = 0; // number of multipolygon and multilinestringrelations
    for (auto ri = rels.begin (); ri != rels.end (); ++ri)
    {

    }
    std::vector <bool> mp_okay (nmp);
    std::fill (mp_okay.begin (), mp_okay.end (), true);

    size_t ncol = unique_vals.k_rel.size ();
    rel_id.reserve (nmp);

    Rcpp::CharacterMatrix kv_mat (Rcpp::Dimension (nmp, ncol));
    unsigned int count = 0;

    for (auto itr = rels.begin (); itr != rels.end (); ++itr)
    {
        Rcpp::checkUserInterrupt ();

        rel_id.push_back (std::to_string (itr->id));
        lon_arr_rel.push_back (lon_vec);
        lat_arr_rel.push_back (lat_vec);
        rowname_arr_mp.push_back (rowname_vec);
        id_vec_rel.push_back (ids_rel);

        if (rowname_vec.size () == 0)
            mp_okay [count] = false;

        lon_vec.clear ();
        lon_vec.shrink_to_fit ();
        lat_vec.clear ();
        lat_vec.shrink_to_fit ();
        rowname_vec.clear ();
        rowname_vec.shrink_to_fit ();
        ids_rel.clear ();
        ids_rel.shrink_to_fit ();

        osm_convert::get_value_mat_rel (itr, unique_vals, kv_mat, count++);
    }

    // Erase any multipolygon ways that are not okay. An example of these is
    // opq("salzburg") %>% add_osm_feature (key = "highway"), for which
    // $osm_multipolygons [[42]] with way#4108738 is not okay.
    std::vector <std::string> not_okay_id;
    for (size_t i = 0; i < mp_okay.size (); i++)
        if (!mp_okay [i])
            not_okay_id.push_back (rel_id [i]);

    for (std::string i: not_okay_id)
    {
        std::vector <std::string>::iterator it =
            std::find (rel_id.begin (), rel_id.end (), i);
        //size_t j = static_cast <size_t> (std::distance (rel_id.begin (), it));
        long int j = std::distance (rel_id.begin (), it);
        lon_arr_rel.erase (lon_arr_rel.begin () + j);
        lat_arr_rel.erase (lat_arr_rel.begin () + j);
        rowname_arr_mp.erase (rowname_arr_mp.begin () + j);
        id_vec_rel.erase (id_vec_rel.begin () + j);
        rel_id.erase (rel_id.begin () + j);

        // Retain static_cast here because there will generally be very few
        // instances of this loop
        size_t st_nrow = static_cast <size_t> (kv_mat.nrow ());
        Rcpp::CharacterMatrix kv_mat_mp2 (Rcpp::Dimension (st_nrow - 1, ncol));
        // k is int for type-compatible Rcpp indexing
        for (int k = 0; k < kv_mat.nrow (); k++)
        {
            if (k < j)
                kv_mat_mp2 (k, Rcpp::_) = kv_mat (k, Rcpp::_);
            else if (k > j)
                kv_mat_mp2 (k - 1, Rcpp::_) = kv_mat (k, Rcpp::_);
        }
        kv_mat = kv_mat_mp2;
    }

    Rcpp::DataFrame kv_df;
    if (rel_id.size () > 0) // only if there are linestrings
    {
        kv_mat.attr ("names") = unique_vals.k_rel;
        kv_mat.attr ("dimnames") = Rcpp::List::create (rel_id, unique_vals.k_rel);
        kv_df = osm_convert::restructure_kv_mat (kv_mat, true);
    } else
        kv_df = R_NilValue;

    // ****** clean up *****
    lon_arr_rel.clear ();
    lon_arr_rel.shrink_to_fit ();
    lat_arr_rel.clear ();
    lat_arr_rel.shrink_to_fit ();
    rowname_arr_mp.clear ();
    rowname_arr_mp.shrink_to_fit ();

    rel_id.clear ();
    rel_id.shrink_to_fit ();

    Rcpp::List ret (1);
    ret [0] = kv_df;
    return ret;
}

//' get_osm_ways
//'
//' Store OSM ways as data.frame objects.
//'
//' @param kv_df Pointer to Rcpp::DataFrame to hold key-value pairs
//' @param way_ids Vector of <osmid_t> IDs of ways to trace
//' @param ways Pointer to all ways in data set
//' @param unique_vals pointer to all unique values (OSM IDs and keys) in data set
//' @param bbox Pointer to the bbox
//'
//' @noRd
void osm_df::get_osm_ways (Rcpp::DataFrame &kv_df,
        const std::set <osmid_t> &way_ids, const Ways &ways,
        const UniqueVals &unique_vals, const Rcpp::NumericVector &bbox)
{
    size_t nrow = way_ids.size (), ncol = unique_vals.k_way.size ();
    std::vector <std::string> waynames;
    waynames.reserve (way_ids.size ());

    Rcpp::CharacterMatrix kv_mat (Rcpp::Dimension (nrow, ncol));
    std::fill (kv_mat.begin (), kv_mat.end (), NA_STRING);
    unsigned int count = 0;
    for (auto wi = way_ids.begin (); wi != way_ids.end (); ++wi)
    {
        //unsigned int count = static_cast <unsigned int> (
        //        std::distance (way_ids.begin (), wi));
        Rcpp::checkUserInterrupt ();
        waynames.push_back (std::to_string (*wi));

        auto wj = ways.find (*wi);
        osm_convert::get_value_mat_way (wj, unique_vals, kv_mat, count);
        count++;
    }

    kv_df = R_NilValue;
    if (way_ids.size () > 0)
    {
        kv_mat.attr ("names") = unique_vals.k_way;
        kv_mat.attr ("dimnames") = Rcpp::List::create (waynames, unique_vals.k_way);
        if (kv_mat.nrow () > 0 && kv_mat.ncol () > 0)
            kv_df = osm_convert::restructure_kv_mat (kv_mat, false);
    }
}

//' get_osm_nodes
//'
//' Store OSM nodes as a data.frame
//'
//' @param ptxy Pointer to Rcpp::List to hold the resultant geometries
//' @param kv_df Pointer to Rcpp::DataFrame to hold key-value pairs
//' @param nodes Pointer to all nodes in data set
//' @param unique_vals pointer to all unique values (OSM IDs and keys) in data set
//' @param bbox Pointer to the bbox
//'
//' @noRd
void osm_df::get_osm_nodes (Rcpp::DataFrame &kv_df, const Nodes &nodes,
        const UniqueVals &unique_vals, const Rcpp::NumericVector &bbox)
{
    size_t nrow = nodes.size (), ncol = unique_vals.k_point.size ();

    Rcpp::CharacterMatrix kv_mat (Rcpp::Dimension (nrow, ncol));
    std::fill (kv_mat.begin (), kv_mat.end (), NA_STRING);

    std::vector <std::string> ptnames;
    ptnames.reserve (nodes.size ());
    unsigned int count = 0;
    for (auto ni = nodes.begin (); ni != nodes.end (); ++ni)
    {
        // std::distance requires a static_cast which copies each instance and
        // slows this down by lots of orders of magnitude
        //unsigned int count = static_cast <unsigned int> (
        //        std::distance (nodes.begin (), ni));
        if (count % 1000 == 0)
            Rcpp::checkUserInterrupt ();

        ptnames.push_back (std::to_string (ni->first));
        for (auto kv_iter = ni->second.key_val.begin ();
                kv_iter != ni->second.key_val.end (); ++kv_iter)
        {
            const std::string &key = kv_iter->first;
            unsigned int ndi = unique_vals.k_point_index.at (key);
            kv_mat (count, ndi) = kv_iter->second;
        }
        count++;
    }
    if (unique_vals.k_point.size () > 0)
    {
        kv_mat.attr ("dimnames") = Rcpp::List::create (ptnames, unique_vals.k_point);
        kv_df = osm_convert::restructure_kv_mat (kv_mat, false);
    } else
        kv_df = R_NilValue;

    ptnames.clear ();
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **            THE FINAL RCPP FUNCTION CALLED BY osmdata_data_frame    **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

//' rcpp_osmdata_df
//'
//' Return OSM data in data.frame format
//'
//' @param st Text contents of an overpass API query
//' @return Rcpp::List objects of OSM data
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_osmdata_df (const std::string& st)
{
#ifdef DUMP_INPUT
    {
        std::ofstream dump ("./osmdata-df.xml");
        if (dump.is_open())
        {
            dump.write (st.c_str(), st.size());
        }
    }
#endif

    XmlData xml (st);

    const std::map <osmid_t, Node>& nodes = xml.nodes ();
    const std::map <osmid_t, OneWay>& ways = xml.ways ();
    const std::vector <Relation>& rels = xml.relations ();
    const UniqueVals& unique_vals = xml.unique_vals ();

    std::vector <double> lons, lats;
    std::set <std::string> keyset; // must be ordered!
    Rcpp::List dimnames (0);
    Rcpp::NumericMatrix nmat (Rcpp::Dimension (0, 0));

    /* --------------------------------------------------------------
     * 1. Set up bbox
     * --------------------------------------------------------------*/

    std::vector <std::string> colnames, rownames;
    colnames.push_back ("lon");
    colnames.push_back ("lat");

    Rcpp::NumericVector bbox = rcpp_get_bbox_sf (xml.x_min (), xml.y_min (),
                                              xml.x_max (), xml.y_max ());

    /* --------------------------------------------------------------
     * 2. Extract OSM Relations
     * --------------------------------------------------------------*/

    Rcpp::List kv_df_rels = osm_df::get_osm_relations (rels, unique_vals, bbox);
    kv_df_rels.attr ("class") = "data.frame";


    /* --------------------------------------------------------------
     * 3. Extract OSM ways
     * --------------------------------------------------------------*/

    Rcpp::DataFrame kv_df_ways;
    const std::set <osmid_t> ways_ids;
    osm_df::get_osm_ways (kv_df_ways, ways_ids, ways, unique_vals, bbox);


    /* --------------------------------------------------------------
     * 3. Extract OSM nodes
     * --------------------------------------------------------------*/

    // NOTE: kv_df_points is actually an Rcpp::CharacterMatrix, and the
    // following line *should* construct the wrapped data.frame version with
    // strings not factors, yet this does not work.
    //Rcpp::DataFrame kv_df_points = Rcpp::DataFrame::create (Rcpp::_["stringsAsFactors"] = false);
    Rcpp::DataFrame kv_df_points;
    osm_df::get_osm_nodes (kv_df_points, nodes, unique_vals, bbox);


    /* --------------------------------------------------------------
     * 5. Collate all data
     * --------------------------------------------------------------*/

    Rcpp::List ret (4);
    ret [0] = bbox;
    ret [1] = kv_df_points;
    ret [2] = kv_df_ways;
    ret [3] = kv_df_rels;

    std::vector <std::string> retnames {"bbox", "points_kv", "ways_kv",
        "relations_kv"};
    ret.attr ("names") = retnames;

    return ret;
}
