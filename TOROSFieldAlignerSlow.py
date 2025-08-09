def photoAllignNew(master_file, comparator_file, alligned_file, pixel_cutoff = 20, snr_threshold = 100):

    import os
    import numpy as np
    from numpy.ma.extras import vstack
    from astropy.table import Table, Column
    """
    This function will take different csvs of one field and match stars from a master field to a comparator field in cases
    where fields are moved by pixel amounts.

    Required Arguments
    master_file [str]: File path to the master csv
    comparator_file [str]: File path to the comparator csv
    alligned_file [str]: File path you want the new csv to go to

    Optional Arguments
    pixel_cutoff [int]: How far a star can be from a star on comparator list before computer ignores that possible guess
    snr_threshold [float]: Signal to noise threshold, filters out stars before running calculations assuming Poisson noise
    """

    ###Check if files exist
    if os.path.exists(master_file) == False:

        print(f"Master file {master_file} does not exist!")
        return

    if os.path.exists(comparator_file) == False:
        print(f"Comparator file {comparator_file} does not exist!")
        return


    ###Read the csv files
    master_tbl = Table.read(master_file, format='ascii', converters={'obsid': str}, delimiter=',')
    comparator_tbl = Table.read(comparator_file, format='ascii', converters = {'obsid': str}, delimiter=',')
    print(master_tbl, "MASTER RAW")
    print(comparator_tbl, "COMPARATOR RAW")


    ###Due to background sub, some pixels have negative flux. Let's fix that.
    master_tbl = master_tbl[master_tbl['flux'] >= 0]
    comparator_tbl = comparator_tbl[comparator_tbl['flux'] >= 0]
    print(master_tbl, "MASTER FIXED")
    print(comparator_tbl, "COMPARATOR FIXED")


    ###Will filter out fluxes that dont meet a required SNR assuming poisson noise.
    master_tbl = master_tbl[master_tbl['flux'] / np.sqrt(master_tbl['flux']) >= snr_threshold]
    comparator_tbl = comparator_tbl[comparator_tbl['flux'] / np.sqrt(comparator_tbl['flux']) >= snr_threshold]
    print(master_tbl, "MASTER FILTERED")
    print(comparator_tbl, "COMPARATOR FILTERED")


    ###Next we order the tables from highest flux to lowest
    master_tbl = master_tbl[np.argsort(master_tbl['flux'])[::-1]]
    comparator_tbl = comparator_tbl[np.argsort(comparator_tbl['flux'])[::-1]]
    print(master_tbl, "MASTER SORTED")
    print(comparator_tbl, "COMPARATOR SORTED")

    master_tbl_copy = master_tbl.copy() #Mutable copy


    ###Every frame is slightly shifted by various pixel amounts. To avoid mismatching, we must account for this.
    shift_iteration = 1
    top30_master = master_tbl[np.argsort(master_tbl['flux'])[::-1][:30]]
    top30_comp = comparator_tbl[np.argsort(comparator_tbl['flux'])[::-1][:30]]
    cumulative_translation = np.array([0.0, 0.0])

    while shift_iteration <= 5:

        print("\n" + "=" * 40)
        print("BEGINNING SHIFT CORRECTION STAGE")
        print("=" * 40 + "\n")

        print("Accounting for shift, iteration: ", shift_iteration)

        #Get coords of top 30 brightest stars from both lists
        master_points = vstack([top30_master['xcentroid'] , top30_master['ycentroid']]).T
        comparator_points = vstack([top30_comp['xcentroid'] , top30_comp['ycentroid']]).T
        matched_coords = []

        for c_id , c_point in enumerate(comparator_points):

            distances = np.sqrt((c_point[0] - master_points[:,0])**2 + (c_point[1] - master_points[:,1])**2) #Distances

            candidate_matches = np.where(distances <= pixel_cutoff)[0] #Create candidates list where distance smaller than a threshold

            if len(candidate_matches) == 0:
                continue #Skip this c_id,c_point if no candidate matches

            #Save index of best mpoint match and append to list
            best_mid = candidate_matches[np.argmin(distances[candidate_matches])]
            matched_coords.append((best_mid, c_id))

            print(f"On iteration {shift_iteration} with {len(candidate_matches)} matches")

        #Now we will apply a affine transform and see what matches

        master_coords = np.array( [[top30_master['xcentroid'][i] , top30_master['ycentroid'][i]] for i,_ in matched_coords]) #Array of mtable matches
        comparator_coords = np.array( [[top30_comp['xcentroid'][j] , top30_comp['ycentroid'][j]] for _,j in matched_coords]) #Array of ctable matches

        translation = np.mean(master_coords - comparator_coords , axis = 0) #We get our translation
        cumulative_translation += translation

        #Now we apply this translation to top30M and top30C to repeat to fine tune our shift
        top30_comp['xcentroid'] += translation[0]
        top30_comp['ycentroid'] += translation[1]

        print(f"On iteration {shift_iteration}, translation is {translation} and cumulative translation is {cumulative_translation}")
        shift_iteration += 1

    #Apply our found translation to the comparator table
    comparator_tbl['xcentroid'] += cumulative_translation[0]
    comparator_tbl['ycentroid'] += cumulative_translation[1]





    ###First, we compare all comparator stars to master stars and get the distances and flux differences
    cstar_to_mstar_match_list = []
    used_master_ids = set()
    for cstar in comparator_tbl:

        best_cstar_match = None
        best_cstar_distance = float('inf')
        best_cstar_deltaF = float('inf')

        for mstar in master_tbl:

            if mstar['id'] in used_master_ids:
                continue #if current mstar has already been matched, skip it (matches must be unique)

            delta_x = mstar['xcentroid'] - cstar['xcentroid'] #X shift
            delta_y = mstar['ycentroid'] - cstar['ycentroid'] #Y shift
            delta_F = abs(mstar['flux'] - cstar['flux']) #Flux shift

            distance = np.sqrt((delta_x ** 2) + (delta_y ** 2))

            if distance > 7:

                continue #If this mstar is completely way too far, ignore it!


            #If this mstar has a closer distance and/or it has same distance as best but nearer flux, update cstar match!
            if distance < best_cstar_distance or (distance == best_cstar_distance and delta_F < best_cstar_deltaF):

                best_cstar_match = mstar
                best_cstar_distance = distance
                best_cstar_deltaF = delta_F


        if best_cstar_match is not None:

            used_master_ids.add(best_cstar_match['id']) #Mark this mstar as used
            print(f"Comparator star {cstar['id']} matches with Master star {best_cstar_match['id']} with distance {best_cstar_distance} and ΔF {best_cstar_deltaF}")

        cstar_to_mstar_match_list.append((cstar , best_cstar_match))

    ###Next, we create a reverse comparison of mstars to cstars to keep the matching symmetric
    mstar_to_cstar_match_list = []
    used_comparator_ids = set()
    for mstar in master_tbl:

        best_mstar_match = None
        best_mstar_distance = float('inf')
        best_mstar_deltaF = float('inf')

        for cstar in comparator_tbl:

            delta_x = mstar['xcentroid'] - cstar['xcentroid']
            delta_y = mstar['ycentroid'] - cstar['ycentroid']
            delta_F = abs(mstar['flux'] - cstar['flux'])

            distance = np.sqrt((delta_x ** 2) + (delta_y ** 2))

            if distance > 7:
                continue

            if distance < best_mstar_distance or (distance == best_mstar_distance and delta_F < best_mstar_deltaF):

                best_mstar_match = cstar
                best_mstar_distance = distance
                best_mstar_deltaF = delta_F

        if best_mstar_match is not None:

            used_comparator_ids.add(best_mstar_match['id'])
            print(f"Master star {mstar['id']} matches with comparator star {best_mstar_match['id']} with distance {best_mstar_distance} and ΔF {best_mstar_deltaF}")

        mstar_to_cstar_match_list.append((mstar , best_mstar_match))

    ###Now we will only allow the stars in that choose the same star both ways (symmetric)
    c_m_dict = {c['id'] : m['id'] for (c , m) in cstar_to_mstar_match_list if m is not None} #Matches cstar ids to matched mstar ids
    m_c_dict = {m['id'] : c['id'] for (m , c) in mstar_to_cstar_match_list if c is not None} #Matches mstar ids to matched cstar ids

    final_id_matches = []
    for (c_id , m_id) in c_m_dict.items():

        if m_id in m_c_dict and m_c_dict[m_id] == c_id:

            final_id_matches.append((c_id, m_id))

    print("Total symmetric matches: ", len(final_id_matches))

    ###Now that we have matched a cstar to an mstar, lets build our astropy table again and save the required data!
    aligned_tbl = Table(names=('alligned_id', 'old_id', 'x_n', 'y_n', 'RA', 'DEC', 'f_n'),
                        dtype=('i4', 'i4', 'f8', 'f8', 'S20', 'S20', 'f8')) #Initializing the table

    for (c_id, m_id) in final_id_matches:
        cstar_info = comparator_tbl[comparator_tbl['id'] == c_id][0]

        aligned_tbl.add_row(
            (m_id,
             c_id,
             float(cstar_info['xcentroid'] - cumulative_translation[0]), #for accurate position tracking on every image, need to remove the
             float(cstar_info['ycentroid'] - cumulative_translation[1]), #shift since its an artificial thing.
             str(cstar_info['RA']),
             str(cstar_info['DEC']),
             str(cstar_info['flux']))
        )

    aligned_tbl['Date'] = [comparator_tbl[comparator_tbl['id'] == c]['Date'][0] for (c, _) in final_id_matches] #Add dates
    aligned_tbl['Time'] = [comparator_tbl[comparator_tbl['id'] == c]['Time'][0] for (c, _) in final_id_matches] #Add times

    aligned_tbl.sort('alligned_id') #Sort by id
    aligned_tbl.write(alligned_file, delimiter=',', format='ascii', overwrite=True) #Write to csv
