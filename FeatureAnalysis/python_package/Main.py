import sys
import Config
import os
import DeadOrganoidClassifier
import ProcessFeatures

if __name__ == "__main__":
    try:
        cmd = sys.argv[1]
    except IndexError:
        print "Usage: %s [COMMAND] <options>" % sys.argv[0]
        print "Run '%s help' for a list of commands" % sys.argv[0]
        sys.exit()

    if cmd == "help":
        print "Usage: %s [COMMAND] <options>" % sys.argv[0]
        print "Possible COMMAND <option> combinations:"
        print "- FULL_WORKFLOW <PLATE_ID>"
        print "- NORMALIZE_WELL <PLATE_ID>"
        print "- CLASSIFY_ORGANOIDS <PLATE_ID>"
        print "- AVERAGE_WELLS <PLATE_ID>"
        print "- BLURRY_ORGANOID_STATS <PLATE_ID>"
        print "- DMSO_NORM <PLATE_ID>"
    elif cmd == "FULL_WORKFLOW":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s NORMALIZE_WELL [plate_id]" % sys.argv[0]
            sys.exit()

        # Calc DMSO average
        cls = ProcessFeatures.get_dmso_average_for_plate(plate=plate)

        # Normalize wells
        all_wells = sorted([
            well for well in
            os.listdir(os.path.join(Config.FEATUREDIR, plate, "wells"))
            if well.startswith(plate)])
        for well in all_wells:
            print well
            try:
                FEATURES = ProcessFeatures.load_organoid_features(wells=[well])
                FEATURES = ProcessFeatures.label_blurry_organoids(
                    **FEATURES)
                FEATURES = ProcessFeatures.get_normalized_organoid_features(
                    **FEATURES)
            except:
                continue

        # Calc plate average
        wellavg = ProcessFeatures.calc_well_summaries(plate=plate)

        # Calc blurry organoid stats
        blurry_cls = ProcessFeatures.create_blurry_organoid_statistics(
            plate=plate)

    elif cmd == "NORMALIZE_WELLS":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s NORMALIZE_WELL [plate_id]" % sys.argv[0]
            sys.exit()
        all_wells = sorted([
            well for well in
            os.listdir(os.path.join(Config.FEATUREDIR, plate, "wells"))
            if well.startswith(plate)])
        for well in all_wells:
            print well
            try:
                FEATURES = ProcessFeatures.load_organoid_features(wells=[well])
                FEATURES = ProcessFeatures.label_blurry_organoids(
                    **FEATURES)
                FEATURES = ProcessFeatures.get_normalized_organoid_features(
                    **FEATURES)
            except:
                continue
    elif cmd == "CLASSIFY_ORGANOIDS":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s CLASSIFY_ORGANOIDS [plate_id]" % sys.argv[0]
            sys.exit()
        cls = DeadOrganoidClassifier.classify_organoids(plate=plate)
    elif cmd == "AVERAGE_WELLS":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s AVERAGE_WELLS [plate_id]" % sys.argv[0]
            sys.exit()
        wellavg = ProcessFeatures.calc_well_summaries(plate=plate)
    elif cmd == "BLURRY_ORGANOID_STATS":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s BLURRY_ORGANOID_STATS [plate_id]" % sys.argv[0]
            sys.exit()
        cls = ProcessFeatures.create_blurry_organoid_statistics(plate=plate)
    elif cmd == "DMSO_NORM":
        try:
            plate = sys.argv[2]
        except IndexError:
            print "Usage: %s DMSO_NORM [plate_id]" % sys.argv[0]
            sys.exit()
        cls = ProcessFeatures.get_dmso_average_for_plate(plate=plate)
    else:
        print "Usage: %s [COMMAND] <options>" % sys.argv[0]
        print "Run '%s help' for a list of commands" % sys.argv[0]
        sys.exit()