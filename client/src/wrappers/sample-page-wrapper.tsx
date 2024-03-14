import { useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { ISample } from "../models/view-samples.model";
import { samplesService } from "../services/sample/sample.service";
import SamplePage from "../components/sample-page/sample-page";

const SamplePageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [sample, setSample] = useState<ISample>(undefined);
  const [acmgRuleNames, setAcmgRuleNames] = useState<string[]>([]);

  const loc = useLocation();
  const sampleId = loc.pathname.split("/sample/")[1];

  useEffect(() => {
    if (isLoading) {
      samplesService?.getSample({ sampleId: sampleId }).then((res) => {
        if (res.isSuccess) {
          setSample(res.sample);
          setAcmgRuleNames(res.acmgRuleNames);
          setIsLoading(false);
        } else {
          //TODO: handle error
        }
      });
    }
  }, [isLoading, sampleId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return (
      <SamplePage
        sample={sample}
        acmgRuleNames={acmgRuleNames}
        sampleService={samplesService}
      />
    );
  }
};

export default SamplePageWrapper;
