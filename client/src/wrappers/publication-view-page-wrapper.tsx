import { useLocation } from "react-router-dom";
import PublicationViewPage from "../components/publication-view-page/publication-view-page";
import { publicationService } from "../services/publication/publication.service";
import Loader from "../atoms/loader/loader";
import React, { useEffect, useState } from "react";
import { IPublicationPreview } from "../models/publication-search.model";
import { IGetPublicationsByVariantIdResponse } from "../services/publication/publication.dto";

const PublicationViewPageWrapper: React.FunctionComponent = () => {
  const [isLoading, setIsLoading] = useState(true);
  const [data, setData] = useState<IPublicationPreview[]>(undefined);

  const loc = useLocation();
  const variantId = loc.pathname.split("/publication-view/")[1];

  useEffect(() => {
    if (isLoading) {
      publicationService
        ?.getPublicationsByVariantId({ variantId: variantId })
        .then((res) => {
          if (res.isSuccess) {
            setData(convertPubDates(res.publications));
            setIsLoading(false);
          } else {
            //TODO: handle error
          }
        });
    }
  }, [isLoading, variantId]);

  if (isLoading) {
    return <Loader />;
  } else {
    return <PublicationViewPage variantId={variantId} publications={data} />;
  }

  function convertPubDates(
    publications: IPublicationPreview[]
  ): IPublicationPreview[] {
    let updatedPublications: IPublicationPreview[] = [];

    if (publications) {
      //format publication date
      publications.forEach((publication) => {
        let pub = publication;
        pub.date = new Date(publication.date);
        updatedPublications = updatedPublications.concat(pub);
      });
    }

    return updatedPublications;
  }
};

export default PublicationViewPageWrapper;
