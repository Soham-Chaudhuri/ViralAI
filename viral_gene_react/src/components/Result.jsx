import React, { useEffect, useState, useRef } from "react";
import Navbar from "./Navbar";
import axios from "axios";
import { useLocation } from "react-router-dom";
import ViralFeaturePlots from "./ViralFeaturePlots";
import * as $3Dmol from "3dmol"; // Import 3Dmol.js

const Result = () => {
  const location = useLocation();
  const [queryResult, setQueryResult] = useState([]);
  const [maxValue, setMaxValue] = useState(0);
  const [maxIndex, setMaxIndex] = useState(0);
  const [viralData, setViralData] = useState({});
  const [pdbData, setPdbData] = useState(null);
  const { searchQuery } = location.state || {};
  const classMap = {
    0: "Chikungunya",
    1: "Influenza A",
    2: "Influenza B",
    3: "Influenza C",
    4: "Zika",
  };
  const [dataFetched, setDataFetched] = useState(false);

  useEffect(() => {
    if (!dataFetched && searchQuery) {
      fetchPredictData();
      fetchViralData();
      fetchStructurePDB();
      setDataFetched(true);
    }
  }, [searchQuery, dataFetched]);

  useEffect(() => {
    setMaxIndex(
      queryResult.reduce((iMax, x, i, arr) => (x > arr[iMax] ? i : iMax), 0)
    );
    setMaxValue(queryResult[maxIndex]);
  }, [queryResult, maxIndex]);

  const fetchStructurePDB = async () => {
    try {
      const response = await axios.post("http://localhost:8000/getStructure", {
        sequence: searchQuery,
      });
      console.log("Protein Structure PDB:", response.data);
      setPdbData(response.data.pdb);
    } catch (error) {
      console.error("Error:", error.message);
    }
  };

  const fetchPredictData = async () => {
    try {
      const response = await axios.post("http://localhost:8000/predict", {
        sequence: searchQuery,
      });
      console.log("Response:", response.data.prediction[0]);
      setQueryResult(response.data.prediction[0]);
    } catch (error) {
      console.error("Error:", error.message);
    }
  };

  const fetchViralData = async () => {
    try {
      const response = await axios.post("http://localhost:8000/getData", {
        sequence: searchQuery,
      });
      console.log("Viral Data:", response.data);
      setViralData(response.data);
    } catch (error) {
      console.error("Error:", error.message);
    }
  };

  const render3DModel = () => {
    const viewer = $3Dmol.createViewer("protein-3d-view", {
      defaultcolors: $3Dmol.rasmolElementColors,
      disableMouse: true,
    });
    viewer.addModel(pdbData, "pdb");
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    viewer.zoomTo();
    viewer.render();
  };

  useEffect(() => {
    if (pdbData) {
      render3DModel();
    }
  }, [pdbData]);

  return (
    <>
      {queryResult && viralData && (
        <div className="main-bg w-full lg:min-h-screen min-h-[130vh] bg-zinc-800 text-white p-10 flex-col space-y-4">
          <Navbar />
          <div className="text-center pt-20 py-6">
            <h1 className="mb-4 text-xl font-extrabold leading-none tracking-tight text-gray-900 md:text-3xl lg:text-4xl dark:text-white">
              {maxValue > 0.8
                ? `Virus is a strain of ${classMap[maxIndex]}`
                : "Unknown Virus Strain"}
            </h1>
            <p className="mb-2 text-lg font-normal text-gray-500 lg:text-xl sm:px-16 xl:px-48 dark:text-gray-400">
              Some of its properties are plotted against the predefined viral
              genomes below.
            </p>
          </div>
          <div className="text-center mx-auto">
            <ViralFeaturePlots viralData={viralData} />
          </div>
          <div className="mt-4 text-center">
            <h1 className="mb-4 text-xl font-extrabold leading-none tracking-tight text-gray-900 md:text-3xl lg:text-4xl dark:text-white w-1/2 mx-auto">
              Protein Structure Visualization
            </h1>
            {pdbData ? (
              <div
                id="protein-3d-view"
                style={{
                  width: "50%",
                  height: "400px",
                  position: "relative",
                  borderRadius: "0.5rem",
                  overflow: "hidden",
                  pointerEvents: "none",
                  margin: "1rem auto",
                }}
              ></div>
            ) : (
              <p className="text-white mb-4">Loading protein structure...</p>
            )}
            <p className="mb-2 w-3/4 mx-auto text-lg font-normal text-gray-500 lg:text-xl  dark:text-gray-400">
              The visualization shows the 3D structure of the biomolecule from
              the PDB file, highlighting key features such as the atomic
              arrangement, secondary structure elements (e.g., alpha-helices,
              beta-sheets), and interactions with ligands or cofactors. This
              representation helps in understanding the molecule's shape,
              function, and potential binding sites.
            </p>
          </div>
        </div>
      )}
    </>
  );
};

export default Result;
