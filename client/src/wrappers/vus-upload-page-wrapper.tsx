import { useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { vusService } from "../services/vus/vus.service";
import { IAcmgRule } from "../models/acmg-rule.model";
import VusUploadPage from "../components/vus-upload-page/vus-upload-page";
import { samplesService } from "../services/sample/sample.service";

const VusUploadPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [acmgRules, setAcmgRules] = useState<IAcmgRule[]>(undefined);

  const loc = useLocation();
  const vusId = loc.pathname.split("/vus/")[1];

  useEffect(() => {
    if (isLoading) {
      vusService.getAllAcmgRules().then((res) => {
        if (res.isSuccess) {
          setAcmgRules(res.acmgRules);
          setIsLoading(false);
        } else {
          //TODO: handle error
        }
      });
    }
  }, [isLoading, vusId]);

  if (isLoading || !acmgRules) {
    return <Loader />;
  } else {
    return (
      <VusUploadPage
        sampleService={samplesService}
        vusService={vusService}
        acmgRules={acmgRules ?? []}
      />
    );
  }
};

export default VusUploadPageWrapper;
