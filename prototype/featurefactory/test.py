import pyrem as pr


pol = pr.signal.polygraph_from_pkl("/data/pyrem/eeg_vs_emg_data/mouse_A_24_hours_am_to_am.txt.pkl")
factory = pr.features.feature_factory.FeatureFactory([
                #pr.features.feature_families.PeriodFeatures(),
                pr.features.feature_families.PowerFeatures(),
                pr.features.feature_families.EntropyFeatures(),
                pr.features.feature_families.NonLinearFeatures(),
                pr.features.feature_families.WaveletsFeaturesDB4(),
                pr.features.feature_families.HjorthFeatures()
        ])
df = factory.make_features_for_epochs(pol,30,1)
