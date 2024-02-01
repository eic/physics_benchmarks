{
  inputs.epic-nix.url = "github:veprbl/epic-nix";

  inputs.nixpkgs.follows = "epic-nix/nixpkgs";

  outputs = { self, epic-nix, nixpkgs }:
    let

      inherit (nixpkgs) lib;
      supportedSystems = [ "x86_64-linux" "x86_64-darwin" ];

    in
    {

      devShells = lib.genAttrs supportedSystems
        (system:
          with import nixpkgs {
            inherit system;
            overlays = [ epic-nix.overlays.default ];
          };
          {
            default =
              mkShell {
                buildInputs =
                  (builtins.attrValues epic-nix.packages.${system})
                  ++ (with epic-nix.packages.${system}; [
                    geant4.data.G4EMLOW
                    geant4.data.G4ENSDFSTATE
                    geant4.data.G4ENSDFSTATE
                    geant4.data.G4PARTICLEXS
                    geant4.data.G4PhotonEvaporation
                  ])
                  ++ [
                    nlohmann_json
                    snakemake
                  ];
                shellHook = ''
                  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${lib.makeLibraryPath [ fmt ]}
                  export S3_ACCESS_KEY="eicS3read"
                  export S3_SECRET_KEY="eicS3read"
                '';
              };
          });

    };
}
