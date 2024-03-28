import { useLocation } from "react-router-dom";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import VusPage from "../components/vus-page/vus-page";
import { IVus } from "../models/view-vus.model";
import { vusService } from "../services/vus/vus.service";

const VusPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [vus, setVus] = useState<IVus>(undefined);

  const loc = useLocation();
  const vusId = loc.pathname.split("/vus/")[1];

  useEffect(() => {
    if (isLoading) {
      vusService.getVus({ vusId: parseInt(vusId) }).then((res) => {
        if (res.isSuccess) {
          setVus(res.vus);
          setIsLoading(false);
        } else {
          //TODO: handle error
        }
      });
    }
  }, [isLoading, vusId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return <VusPage vus={vus} />;
  }
};

export default VusPageWrapper;
