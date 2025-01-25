import React from "react";
import { Scatter } from "react-chartjs-2";
import { Chart as ChartJS, Title, Tooltip, Legend, PointElement, LinearScale } from "chart.js";
ChartJS.register(Title, Tooltip, Legend, PointElement, LinearScale);
const getColor = (virus) => {
  const colors = {
    zika: "rgba(255, 99, 132, 0.2)",
    chikungunya: "rgba(54, 162, 235, 0.2)",
    influenza_a: "rgba(75, 192, 192, 0.2)",
    influenza_b: "rgba(153, 102, 255, 0.2)",
    influenza_c: "rgba(255, 159, 64, 0.2)",
    unknown: "rgba(0, 0, 0, 1)",
  };
  return colors[virus] || "rgba(200, 200, 200, 0.6)";
};

const Graph = ({ viralData, featureX, featureY }) => {
  const data = {
    datasets: Object.keys(viralData).map((virus) => ({
      label: virus,
      data: viralData[virus][featureX].map((xVal, i) => ({
        x: xVal,
        y: viralData[virus][featureY][i],
      })),
      backgroundColor: getColor(virus),
    })),
  };
  const scatterOptions = {
      // responsive: true,
      plugins: {
        title: {
          display: true,
          text: `${featureX} vs ${featureY}`,
          font: { size: 18 },
        },
        legend: {
          display: true,
          position: "top",
        },
      },
      scales: {
        x: {
          title: {
            display: true,
            text: `${featureX} label`,
            font: { size: 14 },
          },
          beginAtZero: true,
        },
        y: {
          title: {
            display: true,
            text: `${featureY} label`,
            font: { size: 14 },
          },
          beginAtZero: true,
        },
      },
    };
  return (
    <>
    <div className="w-full h-[400px]">
      <Scatter data={data} options={scatterOptions} />
    </div>
    </>
  )
}

export default Graph
