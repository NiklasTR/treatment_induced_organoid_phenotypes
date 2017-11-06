import sys
import Config
import os
import LoadFeatures
import BlurryOrganoidClassifier
import NormalizeFeatures
import DeadOrganoidClassifier


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
        print "- NORMALIZE_WELL <PLATE_ID>"
        print "- CLASSIFY_ORGANOIDS <PLATE_ID>"
    elif cmd == "NORMALIZE_WELL":
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
                FEATURES = LoadFeatures.load_organoid_features(wells=[well])
                FEATURES = BlurryOrganoidClassifier.remove_blurry_organoids(
                    **FEATURES)
                FEATURES = NormalizeFeatures.get_normalized_organoid_features(
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
    else:
        print "Usage: %s [COMMAND] <options>" % sys.argv[0]
        print "Run '%s help' for a list of commands" % sys.argv[0]
        sys.exit()