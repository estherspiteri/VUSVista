import React, { useEffect, useState } from "react";
import styles from "./App.module.scss";
import LiteraturePreview from "./components/literature-preview/literature-preview";

type AppProps = {};

const App: React.FunctionComponent<AppProps> = () => {
  const [data, setData] = useState(undefined);

  useEffect(() => {
    fetch("http://127.0.0.1:5000/publications/rs730881106", {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response) => {
        return response.text();
      })
      .then((res) => {
        return JSON.parse(res);
      })
      .then((data) => {
        setData(data.data);
      })
      .catch((error) => console.error("error============:", error));
  }, []);

  return (
    <div>
      {data &&
        Object.values(data?.title).map((title, index) => (
          <LiteraturePreview
            pmid={data.pmid[index]}
            title={data.title[index]}
            date={new Date(data.date[index])}
            authors={data.authors[index]}
            journal={data.journal[index]}
            abstract={data.abstract[index]}
            isSupplementaryMaterialMatch={data.is_sup_mat_match[index]}
          />
        ))}
    </div>
  );
};

export default App;
