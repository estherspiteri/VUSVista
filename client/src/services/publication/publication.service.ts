import {
  IGetPublicationsRequest,
  IGetPublicationsResponse,
} from "./publication.dto";

export class PublicationService {
  async getPublications(
    input: IGetPublicationsRequest
  ): Promise<IGetPublicationsResponse> {
    const result: IGetPublicationsResponse = await fetch(
      `/litvar/publications/${input.rsid}`,
      {
        method: "GET",
        headers: {
          "Content-Type": "application/json;charset=UTF-8",
        },
      }
    )
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const publicationService = new PublicationService();
