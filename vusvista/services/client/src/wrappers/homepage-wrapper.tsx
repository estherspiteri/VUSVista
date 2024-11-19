import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { homepageService } from "../services/homepage/homepage.service";
import { IHomepageResponse } from "../services/homepage/homepage.dto";
import { IHomepageData } from "../models/homepage.model";
import HomePage from "../components/home-page/home-page";

const HomePageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [data, setData] = useState<IHomepageData>(undefined);

  useEffect(() => {
    if (isLoading) {
      homepageService.homepage().then((res) => {
        setData(res.data);
        setIsLoading(false);
      });
    }
  }, [isLoading]);

  if (isLoading) {
    return <Loader />;
  } else {
    return <HomePage data={data} />;
  }
};

export default HomePageWrapper;
