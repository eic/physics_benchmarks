#include <iostream>
#include <string>

#include <ROOT/RDataFrame.hxx>

#include <edm4eic/ReconstructedParticleData.h>

int analyze(std::string file)
{
  // open dataframe
  ROOT::RDataFrame df("events", file, {"GeneratedParticles", "ReconstructedChargedParticles"});

  // count total events
  auto count = df.Count();
  if (count == 0) {
    std::cout << "Error: No events found" << std::endl;
    return -1;
  }

  auto n_tracks = [](const std::vector<edm4eic::ReconstructedParticleData> &p) { return (int) p.size(); };

  auto d = df
  .Define("n_tracks_gen", n_tracks, {"GeneratedParticles"})
  .Define("n_tracks_rec", n_tracks, {"ReconstructedChargedParticles"})
  ;

  auto stats_n_tracks_gen = d.Stats("n_tracks_gen");
  auto stats_n_tracks_rec = d.Stats("n_tracks_rec");
  double mean_num_track_thresh = 0.8;
  if (file.find("135to177deg") != std::string::npos) {
    mean_num_track_thresh = 0.6;
  }
  if (stats_n_tracks_rec->GetMean() < mean_num_track_thresh) {
    std::cout << "Error: too few tracks per events (" << stats_n_tracks_rec->GetMean() << ")" << std::endl;
    stats_n_tracks_gen->Print();
    stats_n_tracks_rec->Print();
    return -1;
  }

  // success
  return 0;
}
