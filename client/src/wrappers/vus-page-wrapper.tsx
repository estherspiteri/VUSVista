import { Navigate, useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import VusPage from "../components/vus-page/vus-page";
import { IVus } from "../models/view-vus.model";
import { vusService } from "../services/vus/vus.service";
import { IAcmgRule } from "../models/acmg-rule.model";
import { samplesService } from "../services/sample/sample.service";

const VusPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [vus, setVus] = useState<IVus>(undefined);
  const [acmgRules, setAcmgRules] = useState<IAcmgRule[]>(undefined);

  const loc = useLocation();
  const vusId = loc.pathname.split("/vus/")[1];

  useEffect(() => {
    if (isLoading) {
      vusService.getVus({ vusId: parseInt(vusId) }).then((res) => {
        if (res.isSuccess) {
          setVus(res.vus);
          setAcmgRules(res.acmgRules);
          setIsLoading(false);
        }
      });
    }
  }, [isLoading, vusId]);

  if (isLoading) {
    return <Loader />;
  } else if (vus === null) {
    return <Navigate to="/view-vus" />;
  } else {
    return (
      <VusPage
        vus={vus}
        acmgRules={acmgRules}
        vusService={vusService}
        sampleService={samplesService}
      />
    );
  }
};

export default VusPageWrapper;
