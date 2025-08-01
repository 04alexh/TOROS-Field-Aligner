#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdlib>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;


struct MatchedPairs //Structure that creates objects representing stars that are "matched"
{
    int comparator_id;
    int master_id;

};


std::vector<std::vector<int>> matchingAlgorithm
(
   const std::vector<int>& c_ids,
   const std::vector<int>& m_ids,
   const std::vector<double>& c_xs,
   const std::vector<double>& c_ys,
   const std::vector<double>& c_fs,
   const std::vector<double>& m_xs,
   const std::vector<double>& m_ys,
   const std::vector<double>& m_fs,
   const double& dist
)
{
    constexpr double eps = 1e-8;

    ////C TO M STAR COMP
    std::vector<MatchedPairs> cstar_to_mstar_match_list;
    std::unordered_set<int> used_Master_IDs;

    for (size_t i = 0 ; i < c_ids.size() ; ++i)
    {

       
        //Load in comparator star information
        int c_id = c_ids[i];
        double c_x = c_xs[i];
        double c_y = c_ys[i];
        double c_f = c_fs[i];

        int best_cstar_id = -1;
        double best_cstar_dist = std::numeric_limits<double>::infinity();
        double best_cstar_deltaF = std::numeric_limits<double>::infinity();

        for (size_t j = 0; j < m_ids.size(); ++j)
        {

            //Check if master id already used
            int m_id = m_ids[j];
            if (used_Master_IDs.find(m_id) != used_Master_IDs.end()) { continue; } 


            //Load in master star information
            double m_x = m_xs[j];
            double m_y = m_ys[j];
            double m_f = m_fs[j];

            //Define delta pos and flux
            double delta_x = m_x - c_x;
            double delta_y = m_y - c_y;
            double delta_F = abs(m_f - c_f);

            //Define distance
            double distance = sqrt(pow(delta_x, 2) + pow(delta_y, 2));

            //Throw out mstars that are way too far
            if (distance > dist) { continue; }

            //If this mstar has a closer distance and /or it has same distance as best but nearer flux, update cstar match!
            if (distance < best_cstar_dist || (abs(distance - best_cstar_dist) < eps) && delta_F < best_cstar_deltaF)
            {

                best_cstar_id = m_id;
                best_cstar_dist = distance;
                best_cstar_deltaF = delta_F;

            }
        }

       
        if (best_cstar_id != -1)
        {

            //Adds used master id to the used list
            //used_Master_IDs.insert(best_cstar_id);

        }

        //Put the c_id m_id pair into a MatchedPairs object, add that object to the match list
        cstar_to_mstar_match_list.push_back(MatchedPairs{ c_id , best_cstar_id });

    }

    ////M to C STAR COMP (for symmetry)
    std::vector<MatchedPairs> mstar_to_cstar_match_list;
    std::unordered_set<int> used_Comparator_IDs;

    for (size_t a = 0; a < m_ids.size(); ++a)
    {

        //Load in master star information
        int m_id = m_ids[a];
        double m_x = m_xs[a];
        double m_y = m_ys[a];
        double m_f = m_fs[a];

        int best_mstar_id = -1;
        double best_mstar_dist = std::numeric_limits<double>::infinity();
        double best_mstar_deltaF = std::numeric_limits<double>::infinity();

        for (size_t b = 0; b < c_ids.size(); ++b)
        {

            //Check if comparator id already used
            int c_id = c_ids[b];
            if (used_Comparator_IDs.find(c_id) != used_Comparator_IDs.end()) { continue; }


            //Load in comparator star information
            double c_x = c_xs[b];
            double c_y = c_ys[b];
            double c_f = c_fs[b];

            //Define delta pos and flux
            double delta_x = m_x - c_x;
            double delta_y = m_y - c_y;
            double delta_F = abs(m_f - c_f);

            //Define distance
            double distance = sqrt(pow(delta_x, 2) + pow(delta_y, 2));

            //Throw out cstars that are way too far
            if (distance > dist) { continue; }

            //If this cstar has a closer distance and /or it has same distance as best but nearer flux, update mstar match!

            if (distance < best_mstar_dist || (abs(distance - best_mstar_dist) < eps && delta_F < best_mstar_deltaF))
            {

                best_mstar_id = c_id;
                best_mstar_dist = distance;
                best_mstar_deltaF = delta_F;

            }
        }



        if (best_mstar_id != -1)
        {

            //Adds used comparator id to the used list
            //used_Comparator_IDs.insert(best_mstar_id);

        }

        //Put the m_id c_id pair into a MatchedPairs object, add that object to the match list
        mstar_to_cstar_match_list.push_back(MatchedPairs{ best_mstar_id , m_id });
    }
    

    std::vector<std::vector<int>> final_id_matches;

    //Now we will create dictionaries for the ids in both lists

    std::unordered_map<int , int>c_m_dict;
    for (const auto& pair : cstar_to_mstar_match_list)
    {

        if (pair.master_id != -1)
        {

            c_m_dict[pair.comparator_id] = pair.master_id;
        }
    }

    std::unordered_map<int, int>m_c_dict;
    for (const auto& pair : mstar_to_cstar_match_list)
    {

        if (pair.comparator_id != -1)
        {

            m_c_dict[pair.master_id] = pair.comparator_id;
        }
    }

    //Now we check through the dictionaries and see which items have the same matches

    /*for (const auto& pair : c_m_dict)
    {

        int cid = pair.first;
        int mid = pair.second;

        if (m_c_dict.find(mid) != m_c_dict.end() && m_c_dict[mid] == cid)
        {

            std::vector<int> star_pair = { cid , mid };
            final_id_matches.push_back(star_pair);
        }
    }*/
    std::unordered_set<int> used_cids;
    std::unordered_set<int> used_mids;


    for (const auto& pair : c_m_dict)
    {

        int cid = pair.first;
        int mid = pair.second;

        if (m_c_dict.count(mid) && m_c_dict[mid] == cid)
        {

            if (used_cids.find(cid) == used_cids.end() && used_mids.find(mid) == used_mids.end())
            {

                std::vector<int> star_pair = { cid , mid };
                final_id_matches.push_back(star_pair);

                used_cids.insert(cid);
                used_mids.insert(mid);

            }

        }

    }
 
    
    return final_id_matches;
}


PYBIND11_MODULE(starmatcher, m) {
    m.def("matchingAlgorithm", &matchingAlgorithm, "Star matching function");
}
