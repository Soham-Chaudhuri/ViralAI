import React, { useEffect, useState } from "react";
import Navbar from "./Navbar";
import defaultGene from "./defaultgene";
import axios from "axios";
import { useNavigate } from "react-router-dom";
const Home = () => {
  const [searchQuery, setSearchQuery] = useState("");
  const navigate = useNavigate();
  const copyToClipboard = () => {
    navigator.clipboard
      .writeText(defaultGene)
      .then(() => alert("Sequence copied to clipboard!"))
      .catch((err) => console.error("Failed to copy: ", err));
  };
  const searchBtn = (e) => {
    e.preventDefault();
    console.log(searchQuery);
    if(searchQuery.length === 0) {
      alert('Please enter a gene sequence');
      return;
    }
    setSearchQuery(searchQuery.replace(/\n/g, '').replace(/\s+/g, ''));
    navigate('/result', { state: { searchQuery } });
  };
  return (
    <>
      <div className="main-bg w-full lg:min-h-screen min-h-[130vh] bg-zinc-800 text-white p-10 flex-col space-y-4">
        <Navbar />
        <div className="text-center pt-20 py-6">
          <h1 className="mb-4 text-4xl font-extrabold leading-none tracking-tight text-gray-900 md:text-5xl lg:text-6xl dark:text-white">
            Gene Analyzer
          </h1>
          <p className="mb-2 text-lg font-normal text-gray-500 lg:text-xl sm:px-16 xl:px-48 dark:text-gray-400">
            Predict unknown gene sequences and analyze their behavior against
            predefined bacterial and viral genomes with cutting-edge accuracy.
          </p>
          <p className="text-lg font-normal text-gray-500 lg:text-xl sm:px-16 xl:px-48 dark:text-gray-400">
            Our tool analyzes gene sequences by plotting properties such as GC
            content, hydrophilicity, sequence entropy, and more. It also
            provides insights into potential protein structures based on these
            properties, comparing them against predefined sequences for a
            comprehensive analysis.
          </p>
          <button
            type="button"
            className="text-white absolute end-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2 dark:bg-blue-600 dark:hover:bg-blue-700 dark:focus:ring-blue-800"
            onClick={copyToClipboard}
          >
            Copy Sample Sequence
          </button>
        </div>
        <form className="max-w-md mx-auto">
          <label
            htmlFor="default-search"
            className="mb-2 text-sm font-medium text-gray-900 sr-only dark:text-white"
          >
            Search
          </label>
          <div className="relative">
            <div className="absolute inset-y-0 start-0 flex items-center ps-3 pointer-events-none">
              <svg
                className="w-4 h-4 text-gray-500 dark:text-gray-400"
                aria-hidden="true"
                xmlns="http://www.w3.org/2000/svg"
                fill="none"
                viewBox="0 0 20 20"
              >
                <path
                  stroke="currentColor"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                  strokeWidth="2"
                  d="m19 19-4-4m0-7A7 7 0 1 1 1 8a7 7 0 0 1 14 0Z"
                />
              </svg>
            </div>
            <input
              type="search"
              id="default-search"
              className="block w-full p-4 ps-10 text-sm text-gray-900 border border-gray-300 rounded-lg bg-gray-50 focus:ring-blue-500 focus:border-blue-500 dark:bg-gray-700 dark:border-gray-600 dark:placeholder-gray-400 dark:text-white dark:focus:ring-blue-500 dark:focus:border-blue-500"
              placeholder="Search Gene Sequences"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              required
            />
            <button
              type="submit"
              className="text-white absolute end-2.5 bottom-2.5 bg-blue-700 hover:bg-blue-800 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-4 py-2 dark:bg-blue-600 dark:hover:bg-blue-700 dark:focus:ring-blue-800"
              onClick={searchBtn}
            >
              Search
            </button>
          </div>
        </form>
      </div>
    </>
  );
};

export default Home;
