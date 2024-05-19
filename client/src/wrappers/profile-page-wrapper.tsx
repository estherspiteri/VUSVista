import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { IProfile } from "../models/profile.model";
import ProfilePage from "../components/profile-page/profile-page";
import { profileService } from "../services/profile/profile.service";

const ProfilePageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [data, setData] = useState<IProfile>(undefined);

  useEffect(() => {
    if (isLoading) {
      profileService.profile().then((res) => {
        setData(res.profile);
        setIsLoading(false);
      });
    }
  }, [isLoading]);

  if (isLoading) {
    return <Loader />;
  } else {
    return <ProfilePage profile={data} />;
  }
};

export default ProfilePageWrapper;
