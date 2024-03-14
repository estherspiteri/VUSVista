import React from "react";
import styles from "./sample-page.module.scss";
import { ISample } from "../../models/view-samples.model";
import SampleInfo from "./sample-info/sample-info";
import { SampleService } from "../../services/sample/sample.service";

type SamplePageProps = {
  sample: ISample;
  acmgRuleNames: string[];
  sampleService: SampleService;
};

const SamplePage: React.FunctionComponent<SamplePageProps> = (
  props: SamplePageProps
) => {
  return (
    <div className={styles["sample-page-container"]}>
      <div className={styles.title}>Sample Information</div>
      <div className={styles.description}>
        <p>Below you can find information about the sample.</p>
      </div>

      <SampleInfo
        sample={props.sample}
        acmgRuleNames={props.acmgRuleNames}
        sampleService={props.sampleService}
      />
    </div>
  );
};

export default SamplePage;
