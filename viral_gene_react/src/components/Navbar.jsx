import { useGSAP } from "@gsap/react";
import gsap from "gsap";
import React, { useEffect, useRef, useState } from "react";
import { RiMenu3Line } from "react-icons/ri";
import { IoCloseSharp } from "react-icons/io5";
const Navbar = () => {
  const [windowWidth, setWindowWidth] = useState(window.innerWidth);
  const [visible, setVisible] = useState(false);
  const gsapRef = useRef();
  const toggleVisibility = () => {
    setVisible(!visible);
  };
  const handleResize = () => {
    setWindowWidth(window.innerWidth);
  };
  useEffect(() => {
    window.addEventListener("resize", handleResize);
    return () => {
      window.removeEventListener("resize", handleResize);
    };
  }, []);
  useGSAP(() => {
    gsap.from(".animated-div li", {
      opacity: 0,
      x: "-100%",
      duration: 0.5,
      stagger:0.3
    });
  }, [visible]);
  return (
    <>
      <nav className="w-4/5 bg-zinc-600 bg-opacity-50 rounded-3xl z-[99] backdrop-blur-[10px] fixed top-3 left-1/2 -translate-x-1/2 ">
        <div className="flex flex-wrap items-center justify-around">
          <div className="p-3 text-4xl cursor-pointer font-bold tracking-wide">
            <a href="#">
              VIRAL<span className="text-sky-400">AI</span>
            </a>
          </div>
          {windowWidth > 550 ? (
            <div className="p-5">
              <ul className="flex flex-wrap list-none gap-10 font-semibold uppercase">
                <li className="text-xl cursor-pointer group hover:text-sky-400">
                  <a href="#">Home</a>
                  <span className="block max-w-0 group-hover:max-w-full transition-all duration-300 h-0.5 bg-sky-400"></span>
                </li>
                <li className="text-xl cursor-pointer group hover:text-sky-400">
                <a href="#" target="_blank" rel="noopener noreferrer">Viruses</a>
                  <span className="block max-w-0 group-hover:max-w-full transition-all duration-300 h-0.5 bg-sky-400"></span>
                </li>
                <li className="text-xl cursor-pointer group hover:text-sky-400">
                <a href="#" target="_blank" rel="noopener noreferrer">GitHub</a>
                  <span className="block max-w-0 group-hover:max-w-full transition-all duration-300 h-0.5 bg-sky-400"></span>
                </li>
                <li className="text-xl cursor-pointer group hover:text-sky-400">
                  <a href="#">About</a>
                  <span className="block max-w-0 group-hover:max-w-full transition-all duration-300 h-0.5 bg-sky-400"></span>
                </li>
              </ul>
            </div>
          ) : (
            <>
            {!visible?(
                <RiMenu3Line size={30} onClick={toggleVisibility} />
            ):(
                <IoCloseSharp size={30} onClick={toggleVisibility} />
            )}
              {visible && (
                // <div >
                // </div>
                  <div ref={gsapRef} className="px-5 pt-2 pb-5 animated-div">
                    <ul className="flex-col list-none font-semibold uppercase text-center">
                      <li className="text-xl cursor-pointer group hover:text-sky-400 py-3">
                        <a href="#">Home</a>
                      </li>
                      <li className="text-xl cursor-pointer group hover:text-sky-400 py-3">
                      <a href="#" target="_blank" rel="noopener noreferrer">Viruses</a>
                      </li>
                      <li className="text-xl cursor-pointer group hover:text-sky-400 py-3">
                      <a href="#" target="_blank" rel="noopener noreferrer">GitHub</a>
                      </li>
                      <li className="text-xl cursor-pointer group hover:text-sky-400 py-3">
                        <a href="#">About</a>
                      </li>
                    </ul>
                  </div>
              )}
            </>
          )}
        </div>
      </nav>
    </>
  );
};

export default Navbar;