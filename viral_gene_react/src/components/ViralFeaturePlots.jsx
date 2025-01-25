import React from "react";
import Graph from "./Graph";
const features = [
  "gc_score",
  "at_score",
  "molecular_weight",
  "hydrophobic_score",
  "hydrophilic_score",
  "sequence_entropy",
];
const properties = ["The AT and GC content are important because they can influence the stability of the DNA molecule. GC pairs form three hydrogen bonds, making them stronger and more stable than AT pairs, which form two hydrogen bonds. Consequently, regions of DNA with high GC content are typically more stable and may have different melting temperatures compared to regions with high AT content.","In proteins, hydrophobic and hydrophilic regions influence the folding and structure. Hydrophobic residues tend to be buried in the protein's interior to minimize exposure to water, whereas hydrophilic residues are more likely to be on the protein's exterior, interacting with the aqueous environment.","The relationship between sequence entropy and molecular weight. Sequence entropy measures the variability in the composition of a sequence, with higher entropy indicating greater diversity at each position. Molecular weight represents the total mass of the sequence, determined by the sum of the atomic weights of its constituent nucleotides or amino acids. The graph highlights how the complexity of a sequence, reflected by its entropy, correlates with its size, providing insights into the interplay between sequence structure and mass."]
const ViralFeaturePlots = ({ viralData }) => {
  const featurePairs = [];
  //   for (let i = 0; i < features.length; i++) {
  //     for (let j = i + 1; j < features.length; j++) {
  //       featurePairs.push([features[i], features[j]]);
  //     }
  //   }
  featurePairs.push(["gc_score", "at_score"]);
  featurePairs.push(["hydrophobic_score", "hydrophilic_score"]);
  featurePairs.push(["sequence_entropy", "molecular_weight"]);
  return (
    <div className="flex flex-col">
      {featurePairs.map(([featureX, featureY], index) => (
        <>
          <div
            key={index}
            className="scatter-plot-container bg-white my-4 rounded-lg px-4 py-2 w-1/2 mx-auto"
          >
            <Graph
              viralData={viralData}
              featureX={featureX}
              featureY={featureY}
            />
          </div>
          <p className="mb-8 w-3/4 mx-auto text-lg font-normal text-gray-500 lg:text-xl  dark:text-gray-400">
            {properties[index]}
          </p>
        </>
      ))}
    </div>
  );
};

export default ViralFeaturePlots;
