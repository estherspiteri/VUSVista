import React, { useEffect, useState } from "react";
import styles from "./view-samples-page.module.scss";
import Loader from "../../atoms/loader/loader";
import SampleTable from "./sample-table/sample-table";
import { ISampleSummary } from "../../models/view-samples.model";
import { SampleService } from "../../services/sample/sample.service";

type ViewAllSamplesProps = { sampleService: SampleService };

const ViewAllSamples: React.FunctionComponent<ViewAllSamplesProps> = (
  props: ViewAllSamplesProps
) => {
  const [sampleList, setSampleList] = useState<ISampleSummary[]>(undefined);
  useEffect(() => {
    props.sampleService?.loadAllSamples().then((res) => {
      if (res.isSuccess) {
        setSampleList(res.sampleList);
      }
    });
  }, [props.sampleService]);

  return (
    <div className={styles["view-all-samples-container"]}>
      <div className={styles.title}>Sample List</div>
      <div className={styles.description}>
        <p>
          Below you can find a list of all the samples stored within our
          database.
        </p>
      </div>
      {sampleList ? (
        <div className={styles["samples-container"]}>
          <SampleTable sampleList={sampleList} />
        </div>
      ) : (
        <Loader />
      )}
    </div>
  );
};

export default ViewAllSamples;
